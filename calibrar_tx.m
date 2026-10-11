function res = calibrar_tx(serial, fc, tx_amp)
%CALIBRAR_TX Barre la ganancia del TRANSMISOR y encuentra el nivel correcto.
%
%   En esta bladeRF el control de ganancia RX no tiene efecto (verificado:
%   una senal analogica fuerte no cambia de nivel al mover rx_gain 72 dB).
%   Por tanto el UNICO mando de nivel es el transmisor. Este script lo barre
%   con las antenas colocadas y mide lo que de verdad importa:
%
%     pico|rx|   nivel de la senal recibida. Objetivo: 0.25 - 0.50.
%                Por debajo de 0.05 el eco esta en el ruido.
%                Por encima de 0.7 hay riesgo de saturar en los picos.
%     SNR bin    relacion tono/ruido MEDIDA POR BIN (ambas magnitudes en la
%                misma unidad; comparar un bin contra el ruido de banda
%                completa da un error de 20*log10(sqrt(N)) = 40 dB).
%
%   Regla para el radar de signos vitales:
%     La respiracion necesita ~15 dB de SNR. El latido produce una
%     modulacion de fase 25-30 dB menor, asi que necesita 40-45 dB.
%     Si no llegas a 35 dB, el BPM no va a enganchar de forma fiable.
%
%   IMPORTANTE: pon las antenas en la posicion REAL de medida (25 cm de
%   separacion, apuntando al torax) y a un sujeto o un reflector a 70 cm.
%   El acoplamiento depende muchisimo de la geometria.
%
%   Uso:
%       res = calibrar_tx('8c14dec0ebc947b88424691d0dbf744f');
%       res = calibrar_tx(serial, 5.94e9, 0.15);

if nargin < 1, error('Uso: calibrar_tx(''<serial>'')'); end
if nargin < 2 || isempty(fc),     fc = 5.94e9;  end
if nargin < 3 || isempty(tx_amp), tx_amp = 0.15; end

fs        = 1e6;
bw        = 1e6;
tono      = 100e3;
block_len = 10000;
n_bloques = 6;
rx_gain   = 30;          % da igual: no tiene efecto. Se fija por higiene.

assert(mod(block_len, fs/tono) == 0, 'block_len debe ser multiplo de fs/tono.');

fprintf('Abriendo bladeRF...  fc = %.3f GHz | tx_amplitude = %.3f\n', fc/1e9, tx_amp);
fprintf('Coloca las antenas en la posicion REAL de medida antes de seguir.\n');

b = bladeRF(['*:serial=' serial]);
limpieza = onCleanup(@() cerrar(b));

b.tx.channel = 'TX1';   b.rx.channel = 'RX2';
b.tx.frequency = fc;    b.rx.frequency = fc;
b.tx.samplerate = fs;   b.rx.samplerate = fs;
b.tx.bandwidth = bw;    b.rx.bandwidth = bw;
try, b.rx.agc = 'manual'; catch, end
try, b.rx.gain = rx_gain; catch, end

n        = (0:block_len-1).';
ref      = exp(1j*2*pi*tono*n/fs);
tx_block = tx_amp * ref;
tx_zeros = complex(zeros(block_len,1));

% Bins: tono, DC/fuga de LO, y el resto (ruido)
k_tono = round(tono/fs*block_len) + 1;
k_win  = (k_tono-10):(k_tono+10);
k_dc   = [1:6, block_len-4:block_len];
mask   = true(block_len,1);
mask(k_win) = false;
mask(k_dc)  = false;

b.rx.start();
b.tx.start();
pause(0.3);
for k = 1:30
    b.transmit(tx_block, 5000, 0, false, false);
    b.receive(block_len);
end

tx_gains = [-24 -12 0 10 20 30 40 50 60];
res = struct('tx_gain',[], 'tx_gain_real',[], 'pico',[], ...
             'tono_dbfs',[], 'ruido_bin_dbfs',[], 'snr_db',[], 'pct_clip',[]);

fprintf('\n=== BARRIDO DE GANANCIA TX ===\n');
fprintf('%8s %8s %9s %11s %12s %9s %9s\n', ...
    'tx_gain', 'real', 'pico|rx|', 'tono dBFS', 'ruido/bin', 'SNR dB', '%clip');
fprintf('%s\n', repmat('-', 1, 76));

for i = 1:numel(tx_gains)
    g = tx_gains(i);
    try
        b.tx.gain = g;
    catch err
        fprintf('%8.0f  ERROR: %s\n', g, err.message);
        continue;
    end
    pause(0.20);
    try, g_real = b.tx.gain; catch, g_real = NaN; end

    for k = 1:8
        b.transmit(tx_block, 5000, 0, false, false);
        b.receive(block_len);
    end

    a_tono = []; a_ruido = []; pico = 0; nclip = 0; ntot = 0;
    for k = 1:n_bloques
        b.transmit(tx_block, 5000, 0, false, false);
        [rx, ~, nc, ~] = b.receive(block_len);
        if nc < block_len, continue; end
        seg = rx(1:block_len);

        S = fft(seg);
        % Ambas magnitudes POR BIN, en la misma unidad
        a_tono(end+1,1)  = max(abs(S(k_win))) / block_len;              %#ok<AGROW>
        a_ruido(end+1,1) = sqrt(mean(abs(S(mask)).^2)) / block_len;     %#ok<AGROW>

        pico  = max(pico, max(abs(seg)));
        nclip = nclip + nnz(abs(real(seg)) >= 0.99 | abs(imag(seg)) >= 0.99);
        ntot  = ntot + block_len;
    end
    if numel(a_tono) < 3, continue; end

    t_db = 20*log10(median(a_tono)  + 1e-15);
    r_db = 20*log10(median(a_ruido) + 1e-15);
    snr  = t_db - r_db;
    pct  = 100*nclip/max(ntot,1);

    marca = '';
    if pct > 0.01,                       marca = '  <-- SATURA'; end
    if isempty(marca) && pico >= 0.25 && pico <= 0.50, marca = '  <-- buen nivel'; end
    fprintf('%8.0f %8.1f %9.4f %11.1f %12.1f %9.1f %8.3f%%%s\n', ...
        g, g_real, pico, t_db, r_db, snr, pct, marca);

    res.tx_gain(end+1,1)        = g;       %#ok<AGROW>
    res.tx_gain_real(end+1,1)   = g_real;  %#ok<AGROW>
    res.pico(end+1,1)           = pico;    %#ok<AGROW>
    res.tono_dbfs(end+1,1)      = t_db;    %#ok<AGROW>
    res.ruido_bin_dbfs(end+1,1) = r_db;    %#ok<AGROW>
    res.snr_db(end+1,1)         = snr;     %#ok<AGROW>
    res.pct_clip(end+1,1)       = pct;     %#ok<AGROW>
end

b.transmit(tx_zeros, 5000, 0, false, true);

% ---------------------------------------------------------------------
fprintf('\n=== VEREDICTO ===\n');

if numel(res.snr_db) < 3
    fprintf(2, 'Datos insuficientes.\n');
else
    span = max(res.snr_db) - min(res.snr_db);
    fprintf('El SNR varia %.1f dB al mover tx_gain %.0f dB.\n', ...
        span, max(res.tx_gain)-min(res.tx_gain));
    if span < 10
        fprintf(2, ['[X] El TX tampoco responde. Revisa que la antena TX este\n' ...
                    '    conectada a TX1 y que el cable no este suelto.\n']);
    else
        fprintf('[OK] La ganancia TX si controla el nivel.\n');
    end

    % Mejor punto: pico en ventana segura, y de esos el de mayor SNR
    ok = res.pct_clip < 0.01 & res.pico <= 0.55;
    if any(ok)
        idx = find(ok);
        [mejor_snr, j] = max(res.snr_db(idx));
        sel = idx(j);
        fprintf('\nRecomendado: cfg.tx_gain = %.0f;   (pico %.3f, SNR %.1f dB)\n', ...
            res.tx_gain(sel), res.pico(sel), mejor_snr);

        if mejor_snr < 35
            fprintf(2, ['\n[!] SNR maximo de %.1f dB. Suficiente para RESPIRACION,\n' ...
                        '    insuficiente para LATIDO (necesita ~40 dB porque su\n' ...
                        '    modulacion de fase es 25-30 dB menor).\n' ...
                        '    Para ganar SNR: acerca el sujeto, mejora el apuntado,\n' ...
                        '    o sube tx_amplitude (ahora %.2f, tope 0.5).\n'], ...
                        mejor_snr, tx_amp);
        else
            fprintf('[OK] SNR de %.1f dB: suficiente para intentar el latido.\n', mejor_snr);
        end
    else
        fprintf(2, 'Todos los puntos saturan o pasan de 0.55. Baja tx_amplitude.\n');
    end
end

figure('Name','Calibracion de ganancia TX');
subplot(2,1,1);
plot(res.tx_gain, res.snr_db, '-o', 'LineWidth', 1.4); hold on;
yline(15, ':', 'minimo respiracion');
yline(40, '--', 'minimo latido');
grid on; xlabel('tx\_gain (dB)'); ylabel('SNR por bin (dB)');
title('SNR del eco frente a la ganancia TX');

subplot(2,1,2);
plot(res.tx_gain, res.pico, '-s', 'LineWidth', 1.4); hold on;
yline(0.50, '--r', 'techo seguro');
yline(0.25, '--g', 'suelo util');
grid on; xlabel('tx\_gain (dB)'); ylabel('pico |rx|');
title('Nivel recibido -- busca la zona entre las dos lineas');
end

function cerrar(b)
try b.transmit(complex(zeros(10000,1)), 5000, 0, false, true); catch, end
try b.tx.stop(); catch, end
pause(0.05);
try b.rx.stop(); catch, end
clear b;
end
