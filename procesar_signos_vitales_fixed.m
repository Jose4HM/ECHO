function [rpm_est, bpm_est] = procesar_signos_vitales_fixed(buffer_IQ, Fs, dHP, dLP_rr)
% PROCESAR_SIGNOS_VITALES  Extrae respiración y ritmo cardíaco de la señal
% IQ de un radar CW Doppler (bladeRF 2.0).
%
% Entradas:
%   buffer_IQ  – Vector complejo (Nx1) a fs_bb (p. ej. 100 Hz)
%   Fs         – Frecuencia de muestreo del buffer (Hz)
%   dHP        – Filtro paso-alto IIR ya diseñado (quita DC lento)
%   dLP_rr     – Filtro paso-bajo IIR ya diseñado (banda de respiración)
%
% Salidas:
%   rpm_est – Tasa respiratoria estimada (rpm)
%   bpm_est – Frecuencia cardíaca estimada (bpm)

pad_mult = 8;   % Zero-padding para interpolación espectral fina

% ================================================================
% 1. ACONDICIONAMIENTO IQ
% ================================================================
I = real(buffer_IQ);
Q = imag(buffer_IQ);

% Quitar DC de cada canal
I = I - mean(I);
Q = Q - mean(Q);

% Corrección de desequilibrio I/Q (Gram-Schmidt)
alpha   = dot(I, Q) / (dot(I, I) + eps);
Q_ortho = Q - alpha * I;
g       = std(I) / (std(Q_ortho) + eps);
Q_ortho = Q_ortho * g;
S_iq    = I + 1j * Q_ortho;

% ================================================================
% 2. CORRECCIÓN CFO Y EXTRACCIÓN DE FASE
% ================================================================
% ✅ FIX BUG 4: Solo un detrend, aplicado ANTES de unwrap para
% evitar doble sustracción que recorta la señal respiratoria lenta.
Nbb = numel(S_iq);
n   = (0:Nbb-1).';

phi_raw = unwrap(angle(S_iq));

% Estimar y sustraer rampa de CFO (drift lineal de frecuencia)
p       = polyfit(n, phi_raw, 1);
phi_cfo = phi_raw - polyval(p, n);   % Quita tendencia lineal (CFO + offset)

% Filtro paso-alto para quitar drift lento residual
phi_clean = filtfilt(dHP, phi_cfo);
x         = phi_clean(:);

% ================================================================
% 3. ESTIMACIÓN DE RESPIRACIÓN (RR)
% ================================================================
x_rr = filtfilt(dLP_rr, x);

% FFT con ventana de Hamming y zero-padding
Nfft_rr  = 2^nextpow2(pad_mult * numel(x_rr));
bin_hz   = Fs / Nfft_rr;                          % Ancho de bin correcto
X_rr     = fft(x_rr .* hamming(numel(x_rr)), Nfft_rr);
K_rr     = floor(Nfft_rr/2) + 1;
f_rr_ax  = (0:K_rr-1) * bin_hz;
P_rr     = abs(X_rr(1:K_rr));

% Banda fisiológica de respiración: 6–48 RPM = 0.10–0.80 Hz
mask_rr  = (f_rr_ax >= 0.10) & (f_rr_ax <= 0.80);
f_rr_est = NaN;

if any(mask_rr)
    idx_band        = find(mask_rr);
    [~, i_rel]      = max(P_rr(idx_band));
    k               = idx_band(i_rel);

    % ✅ FIX BUG 2: paso de bin correcto = Fs/Nfft_rr (= bin_hz)
    if k > 1 && k < length(P_rr)
        al = P_rr(k-1); be = P_rr(k); ga = P_rr(k+1);
        denom = al - 2*be + ga;
        delta_idx = 0;
        if abs(denom) > 1e-12
            delta_idx = 0.5 * (al - ga) / denom;
        end
    else
        delta_idx = 0;
    end
    f_peak = f_rr_ax(k) + delta_idx * bin_hz;   % ← bin_hz, no el doble

    % Lógica de sub-armónicos: solo bajar si el candidato/2 es fisiológico
    f_candidate = f_peak;
    k_half = round(k / 2);
    if k_half >= 2 && k_half < length(P_rr)
        if P_rr(k_half) / (P_rr(k) + eps) > 0.6
            f_candidate = f_rr_ax(k_half);
        end
    end

    % ✅ FIX BUG 5: guardia de mínimo fisiológico (6 RPM = 0.10 Hz)
    f_set = [f_candidate, f_candidate/2];
    f_set = f_set(f_set >= 0.10);          % Descartar < 6 RPM
    if isempty(f_set), f_set = f_candidate; end

    scores = zeros(size(f_set));
    for i = 1:numel(f_set)
        f0 = f_set(i);
        [~, k0] = min(abs(f_rr_ax - f0));
        [~, k1] = min(abs(f_rr_ax - 2*f0));
        if k1 <= numel(P_rr)
            scores(i) = P_rr(k0) + 0.5 * P_rr(k1);
        else
            scores(i) = P_rr(k0);
        end
    end
    [~, best] = max(scores);
    f_rr_est  = f_set(best(1));   % Garantizar escalar
end

% Fallback si no se encontró nada razonable
if isnan(f_rr_est) || f_rr_est < 0.10 || f_rr_est > 0.80
    f_rr_est = 0.25;   % 15 RPM: valor adulto en reposo típico
end
rpm_est = double(60 * f_rr_est(1));   % Garantizar escalar

% ================================================================
% 4. CANCELACIÓN DE RESPIRACIÓN (SIC — Successive Interference Cancellation)
% ================================================================
% ✅ FIX BUG 1: construir modelo sobre x_rr Y sustraer de x_rr → residuo
% coherente. Luego restar de x (sin filtrar) solo si x_rr ≈ x en banda HR.
tbb = (0:numel(x)-1).' / Fs;
H   = [ones(numel(tbb),1), tbb, ...
    sin(2*pi*  f_rr_est *tbb), cos(2*pi*  f_rr_est *tbb), ...
    sin(2*pi*2*f_rr_est *tbb), cos(2*pi*2*f_rr_est *tbb), ...
    sin(2*pi*3*f_rr_est *tbb), cos(2*pi*3*f_rr_est *tbb)];

% Ajustar theta sobre x_rr (la señal en banda de respiración)
theta = H \ x_rr;

% Restar la respiración estimada de la MISMA señal x_rr para obtener
% residuo limpio en la banda baja, luego pasarlo a x para la banda HR.
% Como x_rr es la versión LP de x, la resta en x produce la señal HR.
x_rrfree = x - H * theta;


% ================================================================
% 5. ESTIMACIÓN DE RITMO CARDÍACO (HR)
% ================================================================
margin  = max(0.04, 0.08 * f_rr_est);
fmin_hr = max(0.65, f_rr_est + margin);
fmax_hr = 3.33;   % 200 BPM máximo razonable
if fmin_hr >= fmax_hr - 0.1
    fmin_hr = 0.65;
end




% Butterworth bandpass en banda cardíaca
[b_hr, a_hr] = butter(4, [fmin_hr, fmax_hr] / (Fs/2), 'bandpass');
x_hr         = filtfilt(b_hr, a_hr, x_rrfree);

Nfft_hr  = 2^nextpow2(pad_mult * numel(x_hr));
bin_hr   = Fs / Nfft_hr;
X_hr     = fft(x_hr .* hamming(numel(x_hr)), Nfft_hr);
K_hr     = floor(Nfft_hr/2) + 1;
f_hr_ax  = (0:K_hr-1) * bin_hr;
P_hr     = abs(X_hr(1:K_hr));


% --- INICIO BLOQUE DE DEBUG ---
% x_hr es la señal cardíaca filtrada y lista para la FFT
N_debug = numel(x_hr);
f_debug = (0:N_debug-1) * (Fs / N_debug);
X_debug = abs(fft(x_hr .* hamming(N_debug), N_debug)); 

subplot(2, 2, [3, 4]);
% figure;
plot(f_debug(1:floor(N_debug/2)), X_debug(1:floor(N_debug/2)));
title('Espectro de frecuencia de la señal cardíaca (x_hr)');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;
% --- FIN BLOQUE DE DEBUG ---

% ✅ FIX BUG 3: usar mask_final (con armónicos excluidos) en lugar
% de mask_valid al construir los candidatos
delta      = max(0.05, 0.10 * f_rr_est);
mask_valid = (f_hr_ax >= fmin_hr) & (f_hr_ax <= fmax_hr);

kill = false(size(f_hr_ax));
for h = 1:3
    kill = kill | (abs(f_hr_ax - h * f_rr_est) <= delta);
end
mask_final = mask_valid & ~kill;
if ~any(mask_final)
    mask_final = mask_valid;   % Fallback si todo queda excluido
end

% ← mask_final, no mask_valid
cands_f = f_hr_ax(mask_final);
cands_P = P_hr(mask_final);

if isempty(cands_f)
    bpm_est = 75;   % Fallback
    return;
end

% Penalización gaussiana sobre armónicos de respiración
pen = zeros(size(cands_f));
for h = 1:3
    pen = pen + exp(-((cands_f - h*f_rr_est).^2) / (2 * delta^2));
end
score = cands_P ./ (1 + 3 * pen);

% Guardia contra NaN/Inf en score (división por cero o cands_P corruptos)
score(~isfinite(score)) = 0;

if ~any(score > 0)
    % Sin candidato válido: devolver fallback fisiológico
    bpm_est = 75;
    return;
end

% max() garantiza escalar; ii(1) cubre el caso raro de empate exacto
[~, ii] = max(score);
ii      = ii(1);
bpm_est = double(60 * cands_f(ii));

% Guardia fisiológica final — también garantiza escalar
bpm_est = double(max(30, min(200, bpm_est(1))));
end
