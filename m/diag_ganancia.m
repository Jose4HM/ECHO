function res = diag_ganancia(serial, etiqueta, fc)
%DIAG_GANANCIA Sondea el control de ganancia RX y caracteriza lo que entra.
%
%   Uso:
%       r1 = diag_ganancia(SERIAL, 'con antena');
%       % ...desenrosca la antena RX...
%       r2 = diag_ganancia(SERIAL, 'sin antena');
%
%   v2: se mide con ESTADISTICA ROBUSTA (mediana y percentiles), no con
%   RMS. Cuando hay interferencia impulsiva que satura, el RMS queda
%   dominado por las muestras railadas -- que valen fondo de escala sin
%   importar la ganancia -- y parece que la ganancia "no hace nada" aunque
%   si este funcionando. La MEDIANA ignora las rafagas y revela el piso de
%   ruido real, que es lo que debe escalar con la ganancia.

if nargin < 1, error('Uso: diag_ganancia(''<serial>'')'); end
if nargin < 2 || isempty(etiqueta), etiqueta = 'sin etiqueta'; end
if nargin < 3 || isempty(fc), fc = 5.8e9; end

fs = 1e6; bw = 1e6; canal = 'RX2';
block_len = 10000; n_bloques = 8;

fprintf('Abriendo bladeRF...  [%s]  fc = %.3f GHz\n', etiqueta, fc/1e9);
b = bladeRF(['*:serial=' serial]);
limpieza = onCleanup(@() cerrar(b));

fprintf('\n=== PROPIEDADES DE b.rx ===\n');
try
    props = properties(b.rx);
    fprintf('  %s\n', strjoin(props(:).', ', '));
catch err
    fprintf('  No se pudieron listar: %s\n', err.message);
end

b.rx.channel = canal;
b.rx.frequency = fc;
b.rx.samplerate = fs;
b.rx.bandwidth = bw;

fprintf('\n=== MODO DE GANANCIA ===\n');
try
    b.rx.agc = 'manual';
    fprintf('  b.rx.agc = ''manual'' -> lee: %s\n', string(b.rx.agc));
catch err
    warning('No se pudo fijar b.rx.agc: %s', err.message);
end

b.rx.start();
pause(0.3);
for k = 1:10, b.receive(block_len); end

ganancias = [-15 0 10 20 30 45 60];
res = struct('etiqueta', etiqueta, 'fc', fc, 'ganancias', ganancias, ...
    'pedida', [], 'real', [], 'mediana', [], 'p99', [], 'pico', [], 'pct_clip', []);

fprintf('\n=== BARRIDO DE GANANCIA (TX APAGADO) ===\n');
fprintf('%7s %7s %10s %10s %8s %9s\n', ...
    'pedida', 'real', 'MEDIANA', 'p99', 'pico', '%clip');
fprintf('%s\n', repmat('-', 1, 58));

for i = 1:numel(ganancias)
    g_ped = ganancias(i);
    try
        b.rx.gain = g_ped;
    catch err
        fprintf('%7.1f  ERROR: %s\n', g_ped, err.message);
        continue;
    end
    pause(0.15);
    try, g_real = b.rx.gain; catch, g_real = NaN; end

    for k = 1:5, b.receive(block_len); end

    a_all = [];
    n_clip = 0; n_tot = 0;
    for k = 1:n_bloques
        [rx, ~, nc, ~] = b.receive(block_len);
        if nc > 0
            seg = rx(1:nc);
            a_all = [a_all; abs(seg)]; %#ok<AGROW>
            n_clip = n_clip + nnz(abs(real(seg)) >= 0.99 | abs(imag(seg)) >= 0.99);
            n_tot = n_tot + nc;
        end
    end
    if isempty(a_all), continue; end

    med  = median(a_all);
    p99  = prctile(a_all, 99);
    pico = max(a_all);
    pct  = 100 * n_clip / max(n_tot, 1);

    fprintf('%7.1f %7.1f %10.5f %10.5f %8.4f %8.3f%%\n', ...
        g_ped, g_real, med, p99, pico, pct);

    res.pedida(end+1,1)   = g_ped;   %#ok<AGROW>
    res.real(end+1,1)     = g_real;  %#ok<AGROW>
    res.mediana(end+1,1)  = med;     %#ok<AGROW>
    res.p99(end+1,1)      = p99;     %#ok<AGROW>
    res.pico(end+1,1)     = pico;    %#ok<AGROW>
    res.pct_clip(end+1,1) = pct;     %#ok<AGROW>
end

% ---------------------------------------------------------------------
% Captura de una traza para ver la FORMA de lo que entra
% ---------------------------------------------------------------------
b.rx.gain = 0; pause(0.15);
for k = 1:5, b.receive(block_len); end
[traza, ~, nc, ~] = b.receive(block_len);
traza = traza(1:nc);
res.traza = traza;

% ---------------------------------------------------------------------
% VEREDICTO basado en la MEDIANA
% ---------------------------------------------------------------------
fprintf('\n=== VEREDICTO [%s] ===\n', etiqueta);

if numel(res.mediana) >= 2 && all(res.mediana > 0)
    rango_med = 20*log10(max(res.mediana) / min(res.mediana));
    rango_g   = max(res.real) - min(res.real);
    fprintf('La mediana varia %.1f dB sobre %.0f dB de ganancia.\n', rango_med, rango_g);
    if rango_med > 0.5 * rango_g
        fprintf('[OK] La ganancia SI se esta aplicando (el piso de ruido escala).\n');
    elseif rango_med < 6
        fprintf(2, '[X] Ni siquiera la mediana escala. La ganancia no llega al hardware.\n');
    else
        fprintf(2, '[!] La mediana escala poco (%.1f dB). Puede haber compresion del front-end.\n', rango_med);
    end
end

if ~isempty(res.pct_clip)
    duty = mean(res.pct_clip);
    fprintf('Duty cycle medio de saturacion: %.3f%%\n', duty);
    if duty > 0.001 && duty < 5
        fprintf(2, ['[X] SATURACION IMPULSIVA (rafagas cortas, no continua).\n' ...
            '    Firma tipica de WiFi de 5 GHz o video FPV. Un tono continuo\n' ...
            '    daria un duty cercano al 100%%, no %.2f%%.\n'], duty);
    elseif duty >= 5
        fprintf(2, '[X] Saturacion casi continua (%.1f%%): hay una portadora fuerte, no rafagas.\n', duty);
    else
        fprintf('[OK] Practicamente sin saturacion.\n');
    end
end

% ---------------------------------------------------------------------
% Grafica de la traza
% ---------------------------------------------------------------------
figure('Name', sprintf('diag_ganancia [%s]', etiqueta));
subplot(2,1,1);
t_ms = (0:numel(traza)-1)/fs*1e3;
plot(t_ms, abs(traza), 'LineWidth', 0.7); hold on;
yline(1.0, ':r', 'fondo de escala');
grid on; xlabel('ms'); ylabel('|x|');
title(sprintf('Envolvente RX, ganancia 0 dB [%s] -- rafagas = interferencia', etiqueta));

subplot(2,1,2);
Nf = 2^nextpow2(numel(traza));
X = fftshift(fft(traza .* hamming(numel(traza)), Nf));
f = (-Nf/2:Nf/2-1) * fs/Nf / 1e3;
plot(f, 20*log10(abs(X)+1e-12), 'LineWidth', 0.8);
grid on; xlabel('kHz desde la portadora'); ylabel('dB');
title(sprintf('Espectro en %.3f GHz (ventana de %.0f kHz)', fc/1e9, fs/1e3));

fprintf('\nSiguiente paso: repite con la antena RX DESENROSCADA:\n');
fprintf('   r2 = diag_ganancia(''%s'', ''sin antena'');\n', serial);
fprintf('Si el pico se desploma -> entra por la antena (interferencia real).\n');
fprintf('Si sigue en 1.41       -> el problema es interno, no del aire.\n');
end

function cerrar(b)
try b.rx.stop(); catch, end
clear b;
end
