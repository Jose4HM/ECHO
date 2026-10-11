function res = validar_captura(archivo)
%VALIDAR_CAPTURA Certifica que la captura es correcta, sin depender del
%   procesamiento de signos vitales.
%
%   ENTRADA: una captura de ESCENA ESTATICA. Sin persona, sin nadie
%   moviendose en la sala, nada oscilando delante. Idealmente con un
%   objeto fijo (una caja, una silla) donde iria el sujeto.
%
%   Sobre una escena estatica, la fase medida DEBERIA ser constante. Lo que
%   se mida de mas es el suelo de ruido del sistema, y es lo que hay que
%   comparar contra la senal fisiologica:
%
%       Respiracion : 4-12 mm    -> 1.0-3.0 rad
%       Latido      : 0.25-0.5 mm -> 0.06-0.12 rad   (a lambda = 50.5 mm)
%
%   Uso:
%       res = validar_captura('estatica.mat');

S = load(archivo);
if isfield(S,'rx_continuous'), x = S.rx_continuous(:); else, x = S.rx(:); end
cfg = S.cfg;
fs  = cfg.fs;
N   = numel(x);
lambda = 3e8 / cfg.fc;

fprintf('\n===== VALIDACION DE CAPTURA: %s =====\n', archivo);
fprintf('%.1f s a %.0f kHz | fc = %.3f GHz | lambda = %.1f mm\n\n', ...
    N/fs, fs/1e3, cfg.fc/1e9, 1e3*lambda);

res = struct();

% ---------------------------------------------------------------- 1
fprintf('--- 1. INTEGRIDAD DEL FLUJO ---\n');
n_clip = nnz(abs(real(x)) >= 0.99 | abs(imag(x)) >= 0.99);
pct_clip = 100*n_clip/N;
fprintf('Saturacion      : %.4f%%\n', pct_clip);
if isfield(cfg,'actual_counts')
    n_inc = nnz(cfg.actual_counts ~= cfg.block_len);
    fprintf('Bloques incompl.: %d de %d\n', n_inc, numel(cfg.actual_counts));
    res.bloques_incompletos = n_inc;
end
res.pct_clip = pct_clip;
ok1 = pct_clip < 0.01;
veredicto(ok1, 'flujo integro', 'hay saturacion: baja tx_amplitude');

% ---------------------------------------------------------------- 2
fprintf('\n--- 2. EL TONO ---\n');
Nf = 2^19;
i0 = max(1, floor(N/2) - Nf/2);
seg = x(i0:i0+Nf-1);
X = fft(seg .* hann(Nf));
fax = (0:Nf-1)'/Nf*fs;  fax(fax>fs/2) = fax(fax>fs/2) - fs;
mag = abs(X);

[~, kt] = max(mag);
f_tone = fax(kt);
% Imagen: misma frecuencia, signo opuesto -> mide el desbalance I/Q
[~, ki] = min(abs(fax + f_tone));
img_db = 20*log10(mag(ki)/mag(kt));

banda = abs(fax - f_tone) > 2e3;   % fuera de +-2 kHz del tono
piso  = sqrt(mean(mag(banda).^2));
snr_db = 20*log10(mag(kt)/piso);

fprintf('Frecuencia      : %.1f Hz (nominal %.0f, offset %+.1f Hz)\n', ...
    f_tone, cfg.tx_tone_hz, f_tone - cfg.tx_tone_hz);
fprintf('SNR del tono    : %.1f dB\n', snr_db);
fprintf('Rechazo imagen  : %.1f dB  (desbalance I/Q)\n', img_db);
res.f_tone = f_tone; res.snr_tono_db = snr_db; res.imagen_db = img_db;
veredicto(snr_db > 25, 'tono solido', 'tono debil: sube tx_gain o acerca el blanco');
veredicto(img_db < -25, 'balance I/Q correcto', 'desbalance I/Q alto (>-25 dB)');

% ---------------------------------------------------------------- 3
fprintf('\n--- 3. CONTINUIDAD ENTRE BLOQUES ---\n');
% cw.m transmite por bloques de cfg.block_len. Si la fase del tono no
% empalma en la frontera, aparece un tren de saltos a fs/block_len = 100 Hz
% que al decimar a 100 Hz se pliega sobre DC y contamina justo la banda
% fisiologica. Se compara el salto de fase EN la frontera contra el salto
% medido a mitad de bloque (control): si son iguales, no hay artefacto.
L = cfg.block_len;
M = 200;
bordes = L:L:(N-L);
if numel(bordes) > 200, bordes = bordes(round(linspace(1,numel(bordes),200))); end

salto_borde  = medir_saltos(x, bordes,       M, f_tone, fs);
salto_centro = medir_saltos(x, bordes + L/2, M, f_tone, fs);

fprintf('Salto en frontera de bloque : %.4f rad (RMS)\n', salto_borde);
fprintf('Salto a mitad de bloque     : %.4f rad (RMS)  <- control\n', salto_centro);
razon = salto_borde / max(salto_centro, 1e-9);
fprintf('Razon frontera/control      : %.2f\n', razon);
res.salto_borde_rad = salto_borde; res.salto_centro_rad = salto_centro;
veredicto(razon < 2, 'sin artefacto de bloque', ...
    'HAY discontinuidad en las fronteras de bloque: el TX no empalma');

% ---------------------------------------------------------------- 4
fprintf('\n--- 4. ESTABILIDAD DE AMPLITUD ---\n');
nb = floor(N/L);
env = zeros(nb,1);
for k = 1:nb, env(k) = rms(x((k-1)*L+1 : k*L)); end
cv = std(env)/mean(env);
caidas = nnz(env < 0.3*median(env));
fprintf('Variacion (CV)  : %.1f%%\n', 100*cv);
fprintf('Caidas profundas: %d bloques por debajo del 30%% de la mediana\n', caidas);
res.cv_amplitud = cv; res.caidas = caidas;
veredicto(caidas == 0, 'amplitud estable', ...
    'hay desvanecimientos: el eco cancela la fuga (separa mas las antenas)');

% ---------------------------------------------------------------- 5
fprintf('\n--- 5. RUIDO DE FASE EN LA BANDA FISIOLOGICA ---\n');
fprintf('    (este es EL numero que decide todo)\n');
cfg2 = cfg; cfg2.frontend_fs_out = 100;
[iq, Fsb] = frontend_dsp(x, cfg2);

z = iq(:);
n = (0:numel(z)-1)';
% Quitar offset de portadora (rampa) para no contarlo como ruido
psi = unwrap(angle(z));
p   = polyfit(n, psi, 3);
phi = psi - polyval(p, n);
% Recortar transitorios
nt = round(2*Fsb);
if numel(phi) > 4*nt, phi = phi(nt+1:end-nt); end

bandas = {'Respiracion 0.15-0.50 Hz', [0.15 0.50], 1.0; ...
          'Cardiaca    0.80-2.00 Hz', [0.80 2.00], 0.093};
res.bandas = struct();
for i = 1:size(bandas,1)
    bp = designfilt('bandpassiir','FilterOrder',6, ...
        'HalfPowerFrequency1',bandas{i,2}(1), ...
        'HalfPowerFrequency2',bandas{i,2}(2),'SampleRate',Fsb);
    r = rms(filtfilt(bp, phi));
    d_um = 1e6 * r * lambda/(4*pi);
    senal = bandas{i,3};
    margen = 20*log10(senal / max(r,1e-12));
    fprintf('%s : %.4f rad = %.1f um  |  margen %+.1f dB\n', ...
        bandas{i,1}, r, d_um, margen);
    res.bandas(i).rms_rad = r;
    res.bandas(i).margen_db = margen;
end

fprintf('\n===== VEREDICTO =====\n');
m_rr = res.bandas(1).margen_db;
m_hr = res.bandas(2).margen_db;
if m_hr > 12
    fprintf('[OK] El LATIDO esta %.0f dB sobre el suelo. La captura sirve para\n', m_hr);
    fprintf('     ambos signos. Si el BPM falla, el problema es el PROCESAMIENTO.\n');
elseif m_hr > 0
    fprintf(2,'[!] El latido queda solo %.0f dB sobre el suelo. Marginal: saldra\n', m_hr);
    fprintf(2,'    a veces. Necesitas mas margen antes de exigirle precision.\n');
else
    fprintf(2,'[X] El latido esta %.0f dB POR DEBAJO del suelo de ruido de fase.\n', -m_hr);
    fprintf(2,'    Ningun algoritmo puede recuperarlo. El limite es la captura,\n');
    fprintf(2,'    no el procesamiento. Deja de tocar el DSP.\n');
end
if m_rr > 12
    fprintf('[OK] La RESPIRACION esta %.0f dB sobre el suelo: solida.\n', m_rr);
else
    fprintf(2,'[!] Incluso la respiracion esta justa (%.0f dB). Revisa el montaje.\n', m_rr);
end
end

% =========================================================================
function s = medir_saltos(x, idx, M, f_tone, fs)
d = zeros(numel(idx),1);
for k = 1:numel(idx)
    b = idx(k);
    ia = b-M+1; ib = b+M;
    if ia < 1 || ib > numel(x), d(k) = NaN; continue; end
    na = (ia-1:b-1)';   nb = (b:ib-1)';
    pa = angle(mean(x(ia:b)    .* exp(-1j*2*pi*f_tone*na/fs)));
    pb = angle(mean(x(b+1:ib)  .* exp(-1j*2*pi*f_tone*nb/fs)));
    d(k) = angle(exp(1j*(pb-pa)));
end
s = rms(d(~isnan(d)));
end

% =========================================================================
function veredicto(ok, msg_ok, msg_mal)
if ok, fprintf('  [OK] %s\n', msg_ok); else, fprintf(2,'  [X] %s\n', msg_mal); end
end
