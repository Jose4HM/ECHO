function tabla = scan_banda(serial, f_ini, f_fin, paso, fs_scan)
%SCAN_BANDA Barre la banda buscando la frecuencia mas limpia para el radar.
%
%   Tus helicoidales cubren 4900-6000 MHz, asi que tienes 1.1 GHz donde
%   elegir. Este script mide cuanta energia ambiental hay en cada tramo,
%   con el TX apagado, y te dice cuales son los mas silenciosos.
%
%   Uso:
%       t = scan_banda('8c14dec0ebc947b88424691d0dbf744f');
%       t = scan_banda(SERIAL, 4.95e9, 5.95e9, 10e6, 10e6);
%
%   Columnas de la tabla:
%       fc_GHz    frecuencia central del tramo
%       mediana   piso de ruido robusto (lo que quieres MINIMO)
%       p99       percentil 99 (sensible a rafagas)
%       pico      maximo observado
%       pct_clip  %% de muestras saturadas (lo que quieres en CERO)
%
%   IMPORTANTE: 5725-5850 MHz es U-NII-3, la banda WiFi mas cargada.
%   5900-6000 MHz suele estar mucho mas limpia y tu antena la cubre.

if nargin < 1, error('Uso: scan_banda(''<serial>'')'); end
if nargin < 2 || isempty(f_ini),   f_ini   = 4.95e9; end
if nargin < 3 || isempty(f_fin),   f_fin   = 5.98e9; end
if nargin < 4 || isempty(paso),    paso    = 10e6;   end
if nargin < 5 || isempty(fs_scan), fs_scan = 10e6;   end

canal = 'RX2';
block_len = 20000;
n_bloques = 3;
ganancia = 0;

frecs = f_ini:paso:f_fin;
n = numel(frecs);

fprintf('Abriendo bladeRF...\n');
b = bladeRF(['*:serial=' serial]);
limpieza = onCleanup(@() cerrar(b));

b.rx.channel = canal;
b.rx.samplerate = fs_scan;
b.rx.bandwidth = min(fs_scan, 20e6);
try, b.rx.agc = 'manual'; catch, end
b.rx.gain = ganancia;
b.rx.frequency = frecs(1);

b.rx.start();
pause(0.3);
for k = 1:10, b.receive(block_len); end

med = zeros(n,1); p99 = zeros(n,1); pico = zeros(n,1); pct = zeros(n,1);

fprintf('\nBarriendo %d tramos de %.1f MHz entre %.2f y %.2f GHz...\n', ...
    n, paso/1e6, f_ini/1e9, f_fin/1e9);
fprintf('%9s %10s %10s %9s %9s\n', 'fc (GHz)', 'mediana', 'p99', 'pico', '%clip');
fprintf('%s\n', repmat('-', 1, 52));

for i = 1:n
    b.rx.frequency = frecs(i);
    pause(0.05);
    for k = 1:4, b.receive(block_len); end   % vaciar tras retunear

    a = []; nclip = 0; ntot = 0;
    for k = 1:n_bloques
        [rx, ~, nc, ~] = b.receive(block_len);
        if nc > 0
            seg = rx(1:nc);
            a = [a; abs(seg)]; %#ok<AGROW>
            nclip = nclip + nnz(abs(real(seg)) >= 0.99 | abs(imag(seg)) >= 0.99);
            ntot = ntot + nc;
        end
    end
    if isempty(a), continue; end

    med(i)  = median(a);
    p99(i)  = prctile(a, 99);
    pico(i) = max(a);
    pct(i)  = 100*nclip/max(ntot,1);

    marca = '';
    if pct(i) > 0.001, marca = '  <-- saturando'; end
    fprintf('%9.3f %10.5f %10.5f %9.4f %8.3f%%%s\n', ...
        frecs(i)/1e9, med(i), p99(i), pico(i), pct(i), marca);
end

tabla = table(frecs(:)/1e9, med, p99, pico, pct, ...
    'VariableNames', {'fc_GHz','mediana','p99','pico','pct_clip'});

% ---- Ranking: primero sin saturacion, luego menor p99 -------------------
limpio = pct < 1e-6;
fprintf('\n=== TRAMOS MAS LIMPIOS ===\n');
if any(limpio)
    sub = tabla(limpio, :);
    [~, ord] = sort(sub.p99, 'ascend');
    sub = sub(ord, :);
    k = min(8, height(sub));
    disp(sub(1:k, :));
    fprintf('Recomendado: cfg.fc = %.4fe9;   %% %.3f GHz\n', ...
        sub.fc_GHz(1), sub.fc_GHz(1));
else
    fprintf(2, 'Ningun tramo quedo libre de saturacion. Reduce la ganancia o\n');
    fprintf(2, 'aleja las antenas de la fuente antes de repetir.\n');
    [~, ord] = sort(tabla.pct_clip, 'ascend');
    disp(tabla(ord(1:min(8,n)), :));
end

% ---- Grafica -------------------------------------------------------------
figure('Name','Barrido de banda 4.9-6.0 GHz');
subplot(2,1,1);
plot(tabla.fc_GHz, 20*log10(tabla.p99 + 1e-12), '-o', 'LineWidth', 1.2); hold on;
plot(tabla.fc_GHz, 20*log10(tabla.mediana + 1e-12), '-s', 'LineWidth', 1.0);
grid on; xlabel('GHz'); ylabel('dBFS');
legend('p99 (rafagas)','mediana (piso)','Location','best');
title('Energia ambiental por tramo (TX apagado)');

subplot(2,1,2);
bar(tabla.fc_GHz, tabla.pct_clip);
grid on; xlabel('GHz'); ylabel('% muestras saturadas');
title('Saturacion por tramo -- busca los tramos en cero');
end

function cerrar(b)
try b.rx.stop(); catch, end
clear b;
end
