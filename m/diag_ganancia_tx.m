function res = diag_ganancia_tx(serial, fc, tx_gain, tx_amp)
%DIAG_GANANCIA_TX Comprueba si la ganancia RX se aplica de verdad, usando el
%   propio transmisor como fuente de referencia conocida.
%
%   v2 - correcciones sobre v1:
%     * El tono se mide por FFT (magnitud), no por correlacion coherente.
%       La v1 promediaba seg.*conj(ref) sobre 10000 muestras: si el buffer
%       de TX sufre un underrun, la fase del tono salta a mitad de bloque y
%       la media se cancela -> lecturas fantasma de -80 dBFS. La FFT no
%       depende de la fase absoluta.
%     * Se ceban DOS bloques de TX por cada recepcion, para que el DAC no
%       se quede sin datos.
%     * Mediana entre bloques en vez de media (un bloque malo no arrastra).
%     * Se reporta la DISPERSION entre bloques. Si es alta, la medida no es
%       fiable y no se emite veredicto.
%     * TX por defecto a un nivel usable. Con tx_gain=-24 y amp=0.03 el
%       tono llega al ruido y todo lo que sigue es azar.
%
%   Uso:
%       diag_ganancia_tx('8c14dec0ebc947b88424691d0dbf744f')
%       diag_ganancia_tx(serial, 5.94e9, 0, 0.10)
%
%   Interpretacion:
%     * El tono sube ~1 dB/dB  -> la ganancia RX funciona.
%     * El tono se queda plano -> el hardware esta clavado en ganancia fija;
%       controla el nivel solo con tx_amplitude / tx_gain.
%     * Dispersion alta       -> medida no valida, no concluyas nada.

if nargin < 1, error('Uso: diag_ganancia_tx(''<serial>'')'); end
if nargin < 2 || isempty(fc),      fc = 5.94e9;  end
if nargin < 3 || isempty(tx_gain), tx_gain = 0;  end   % antes -24 (minimo)
if nargin < 4 || isempty(tx_amp),  tx_amp = 0.10; end  % antes 0.03

fs        = 1e6;
bw        = 1e6;
tono      = 100e3;
block_len = 10000;      % tono cae exactamente en el bin 1000
n_bloques = 8;

assert(mod(block_len, fs/tono) == 0, 'block_len debe ser multiplo de fs/tono.');

fprintf('Abriendo bladeRF...  fc = %.3f GHz | tx_gain = %.0f dB | tx_amp = %.3f\n', ...
    fc/1e9, tx_gain, tx_amp);
b = bladeRF(['*:serial=' serial]);
limpieza = onCleanup(@() cerrar(b));

b.tx.channel = 'TX1';   b.rx.channel = 'RX2';
b.tx.frequency = fc;    b.rx.frequency = fc;
b.tx.samplerate = fs;   b.rx.samplerate = fs;
b.tx.bandwidth = bw;    b.rx.bandwidth = bw;
b.tx.gain = tx_gain;
try, b.rx.agc = 'manual'; catch, end

n         = (0:block_len-1).';
ref       = exp(1j*2*pi*tono*n/fs);
tx_block  = tx_amp * ref;
tx_zeros  = complex(zeros(block_len,1));

% Indices de la FFT
k_tono  = round(tono/fs*block_len) + 1;          % bin del tono (MATLAB 1-based)
k_win   = (k_tono-10):(k_tono+10);               % ventana por si hay deriva
k_dc    = [1:6, block_len-4:block_len];          % DC y fuga de LO
mask    = true(block_len,1);
mask(k_win) = false;
mask(k_dc)  = false;

b.rx.start();
b.tx.start();
pause(0.3);
for k = 1:40
    b.transmit(tx_block, 5000, 0, false, false);
    b.transmit(tx_block, 5000, 0, false, false);
    b.receive(block_len);
end

ganancias = [-15 0 10 20 30 45 60];
res = struct('ganancia_real', [], 'tono_dbfs', [], 'piso_dbfs', [], ...
             'pct_clip', [], 'dispersion_db', []);

fprintf('\n=== NIVEL DEL TONO vs GANANCIA RX (TX ENCENDIDO) ===\n');
fprintf('tono = %.0f kHz, medido por FFT en el bin %d\n\n', tono/1e3, k_tono-1);
fprintf('%7s %7s %11s %11s %9s %10s %9s\n', ...
    'pedida', 'real', 'TONO dBFS', 'piso dBFS', 'SNR dB', 'dispers.', '%clip');
fprintf('%s\n', repmat('-', 1, 72));

for i = 1:numel(ganancias)
    g_ped = ganancias(i);
    try
        b.rx.gain = g_ped;
    catch err
        fprintf('%7.1f  ERROR: %s\n', g_ped, err.message);
        continue;
    end
    pause(0.20);
    try, g_real = b.rx.gain; catch, g_real = NaN; end

    % Purga: descarta lo que ya estuviera en vuelo con la ganancia anterior
    for k = 1:8
        b.transmit(tx_block, 5000, 0, false, false);
        b.transmit(tx_block, 5000, 0, false, false);
        b.receive(block_len);
    end

    a_tono = []; a_piso = []; nclip = 0; ntot = 0;
    for k = 1:n_bloques
        b.transmit(tx_block, 5000, 0, false, false);
        b.transmit(tx_block, 5000, 0, false, false);
        [rx, ~, nc, ~] = b.receive(block_len);
        if nc < block_len, continue; end
        seg = rx(1:block_len);

        S    = fft(seg);
        magA = abs(S) / block_len;                    % amplitud por bin

        a_tono(end+1,1) = max(magA(k_win));           %#ok<AGROW>
        a_piso(end+1,1) = sqrt(sum(abs(S(mask)).^2)) / block_len;  %#ok<AGROW>

        nclip = nclip + nnz(abs(real(seg)) >= 0.99 | abs(imag(seg)) >= 0.99);
        ntot  = ntot + block_len;
    end
    if numel(a_tono) < 3, continue; end

    tono_db = 20*log10(median(a_tono) + 1e-12);
    piso_db = 20*log10(median(a_piso) + 1e-12);
    disp_db = 20*log10(max(a_tono)/max(min(a_tono), 1e-12));
    pct     = 100*nclip/max(ntot,1);

    marca = '';
    if pct > 0.01,    marca = '  <-- satura'; end
    if disp_db > 6,   marca = [marca '  <-- inestable']; end
    fprintf('%7.1f %7.1f %11.1f %11.1f %9.1f %9.1f %8.3f%%%s\n', ...
        g_ped, g_real, tono_db, piso_db, tono_db - piso_db, disp_db, pct, marca);

    res.ganancia_real(end+1,1) = g_real;   %#ok<AGROW>
    res.tono_dbfs(end+1,1)     = tono_db;  %#ok<AGROW>
    res.piso_dbfs(end+1,1)     = piso_db;  %#ok<AGROW>
    res.pct_clip(end+1,1)      = pct;      %#ok<AGROW>
    res.dispersion_db(end+1,1) = disp_db;  %#ok<AGROW>
end

b.transmit(tx_zeros, 5000, 0, false, true);

% ---------------------------------------------------------------------
fprintf('\n=== VEREDICTO ===\n');

ok = res.pct_clip < 0.01 & res.dispersion_db < 6;
n_inest = nnz(res.dispersion_db >= 6);

if n_inest > 0
    fprintf(2, '[!] %d punto(s) con dispersion > 6 dB entre bloques.\n', n_inest);
    fprintf(2, ['    Eso es underrun de TX o interferencia, no ganancia.\n' ...
                '    Se excluyen del ajuste.\n']);
end

snr_med = median(res.tono_dbfs - res.piso_dbfs);
if snr_med < 15
    fprintf(2, ['[X] SNR mediano de solo %.1f dB. El tono esta demasiado cerca\n' ...
                '    del ruido para medir nada. Repite con mas potencia:\n' ...
                '        diag_ganancia_tx(serial, %.2fe9, 10, 0.20)\n'], ...
                snr_med, fc/1e9);
    fprintf(2, '    No interpretes la pendiente hasta tener SNR > 20 dB.\n');
elseif nnz(ok) < 4
    fprintf(2, 'Puntos validos insuficientes (%d). Ajusta tx_gain/tx_amp.\n', nnz(ok));
else
    g = res.ganancia_real(ok);
    t = res.tono_dbfs(ok);
    pend = polyfit(g, t, 1);
    fprintf('Puntos validos           : %d de %d\n', nnz(ok), numel(ok));
    fprintf('Rango de ganancia valido : %.0f dB\n', max(g)-min(g));
    fprintf('Variacion del tono       : %.1f dB\n', max(t)-min(t));
    fprintf('Pendiente                : %.2f dB de tono por dB de ganancia\n', pend(1));
    if pend(1) > 0.7
        fprintf('[OK] La ganancia RX FUNCIONA. Usala normalmente.\n');
    elseif pend(1) < 0.2
        fprintf(2, ['[X] La ganancia RX no tiene efecto: hardware clavado en un\n' ...
                    '    valor fijo. No es bloqueante; controla el nivel solo con\n' ...
                    '    cfg.tx_amplitude y cfg.tx_gain.\n']);
    else
        fprintf(2, '[!] Respuesta parcial (%.2f dB/dB). Posible compresion.\n', pend(1));
    end
end

fprintf(['\nNota: para el radar no elijas rx_gain por este test. Elige la mayor\n' ...
         'ganancia que, con el TX a potencia real y las antenas colocadas,\n' ...
         'deje el pico de |rx| por debajo de 0.5. Eso lo mide cw.m.\n']);

figure('Name','Ganancia RX con tono de referencia (v2)');
plot(res.ganancia_real, res.tono_dbfs, '-o', 'LineWidth', 1.4); hold on;
plot(res.ganancia_real, res.piso_dbfs, '-s', 'LineWidth', 1.1);
g0 = res.ganancia_real(1); t0 = res.tono_dbfs(1);
plot(res.ganancia_real, t0 + (res.ganancia_real - g0), '--k', 'LineWidth', 1.0);
mal = res.dispersion_db >= 6;
if any(mal)
    plot(res.ganancia_real(mal), res.tono_dbfs(mal), 'rx', ...
        'MarkerSize', 12, 'LineWidth', 2);
    legend('tono medido','piso de ruido','ideal 1 dB/dB','medida inestable', ...
        'Location','best');
else
    legend('tono medido','piso de ruido','ideal 1 dB/dB','Location','best');
end
grid on; xlabel('ganancia RX pedida (dB)'); ylabel('dBFS');
title('Si el tono no sigue la linea discontinua, la ganancia no se aplica');
end

function cerrar(b)
try b.transmit(complex(zeros(10000,1)), 5000, 0, false, true); catch, end
try b.tx.stop(); catch, end
pause(0.05);
try b.rx.stop(); catch, end
clear b;
end
