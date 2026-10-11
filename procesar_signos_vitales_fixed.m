function [rpm_est, bpm_est, info] = procesar_signos_vitales_fixed(buffer_IQ, Fs, dHP, dLP_rr)
%PROCESAR_SIGNOS_VITALES_FIXED Estima respiracion y frecuencia cardiaca.
%   Entrada: IQ complejo en banda base (~100 Hz), conservando clutter/DC.
%   La cadena mantiene las correcciones necesarias para bladeRF (offset
%   residual TX/RX y centrado I/Q), pero elimina ramas experimentales que
%   no se usan en la medicion normal.

% ---------------- Parametros ---------------------------------------------
pad_mult    = 8;
trim_s      = 5;
circ_win_s  = 15;
hr_win_s    = 20;
welch_win_s = 15;

rr_min_hz = 6/60;
rr_max_hz = 45/60;
hr_min_hz = 55/60;
hr_max_hz = 160/60;

% ---------------- Validacion ---------------------------------------------
if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('Fs debe ser un escalar positivo en Hz.');
end

buffer_IQ = buffer_IQ(:);
if isempty(buffer_IQ)
    error('buffer_IQ esta vacio.');
end

if any(~isfinite(buffer_IQ))
    idx_good = isfinite(buffer_IQ);
    if nnz(idx_good) < 2
        error('buffer_IQ no tiene suficientes muestras finitas.');
    end
    idx_all = (1:numel(buffer_IQ)).';
    buffer_IQ = interp1(idx_all(idx_good), buffer_IQ(idx_good), ...
        idx_all, 'linear', 'extrap');
end

check_filter_fs(dHP, 'dHP', Fs);
check_filter_fs(dLP_rr, 'dLP_rr', Fs);

dur_cap_s = numel(buffer_IQ)/Fs;
if dur_cap_s < 15
    warning('Captura de %.1f s: RR/HR pueden ser poco estables.', dur_cap_s);
end

if dur_cap_s < 45
    trim_s      = max(1.2, 0.10 * dur_cap_s);
    circ_win_s  = max(5.0, 0.45 * dur_cap_s);
    welch_win_s = min(14.0, max(7.0, 0.40 * dur_cap_s));
    hr_win_s    = max(5.0, 0.40 * dur_cap_s);
end

info = struct( ...
    'rr_fallback', false, ...
    'rr_low_confidence', false, ...
    'rr_prominence', NaN, ...
    'rr_score_ratio', NaN, ...
    'hr_fallback', false, ...
    'hr_low_confidence', false, ...
    'hr_prominence', NaN, ...
    'hr_score_ratio', NaN, ...
    'center_offset_ratio', NaN, ...
    'arc_span_rad', NaN, ...
    'circle_radius', NaN, ...
    'carrier_offset_hz', NaN, ...
    'lo_drift_pp_rad', NaN, ...
    'hr_window_bpm', [], ...
    'hr_agreement', NaN, ...
    'Fs', Fs, ...
    'phi', [], 'x', [], 'x_rr', [], 'x_hr', [], ...
    'f_rr_axis', [], 'P_rr', [], 'f_hr_axis', [], 'P_hr', [], ...
    'f_rr_est', NaN, 'f_hr_est', NaN);

% ========================================================================
% 1. CORRECCION DEL OFFSET RESIDUAL TX/RX
% ========================================================================
S_iq = buffer_IQ;
if dur_cap_s >= 30
    [S_iq, f_off, vueltas, deriva_pp] = quitar_offset_portadora(S_iq, Fs);
else
    [S_iq, f_off, vueltas, deriva_pp] = quitar_offset_portadora(S_iq, Fs, 0);
end
info.carrier_offset_hz = f_off;
info.lo_drift_pp_rad = deriva_pp;

fprintf('\n--- PORTADORA RESIDUAL ---\n');
fprintf('Offset TX-RX estimado       : %+.4f Hz\n', f_off);
fprintf('Vueltas acumuladas          : %.2f en %.1f s\n', vueltas, dur_cap_s);
fprintf('Deriva lenta retirada       : %.2f rad pico a pico\n', deriva_pp);

% ========================================================================
% 2. CENTRADO I/Q Y EXTRACCION DE FASE
% ========================================================================
[c_t, c0, R0, n_fits] = estimar_centro_circulo(S_iq, Fs, circ_win_s);
info.circle_radius = R0;
info.center_offset_ratio = abs(c0)/(R0 + eps);

Z = S_iq - c_t;
phi_raw = unwrap(angle(Z));
info.arc_span_rad = max(phi_raw) - min(phi_raw);

fprintf('\n--- GEOMETRIA IQ ---\n');
fprintf('Radio del arco R            : %.5f\n', R0);
fprintf('|centro| / R                : %.1f\n', info.center_offset_ratio);
fprintf('Apertura angular del arco   : %.2f rad\n', info.arc_span_rad);
fprintf('Ventanas validas de circulo : %d\n', n_fits);

if info.center_offset_ratio > 100
    warning(['|c|/R = %.0f: clutter/fuga domina la senal util. Revise ' ...
        'aislamiento, geometria y orientacion de antenas.'], info.center_offset_ratio);
end
if info.arc_span_rad < 0.25
    warning('Apertura IQ %.2f rad: modulacion toracica muy pequena.', info.arc_span_rad);
end
if info.arc_span_rad > 6
    warning('Apertura IQ %.2f rad: probable movimiento corporal grande.', info.arc_span_rad);
end

% Evitar puntos extremadamente proximos al origen tras el centrado.
n = (0:numel(phi_raw)-1).';
radio = abs(Z);
radio_min = 0.20 * median(radio);
idx_ok = isfinite(phi_raw) & radio > radio_min;
if nnz(idx_ok) >= max(16, round(0.5*numel(phi_raw))) && ~all(idx_ok)
    phi_raw = interp1(n(idx_ok), phi_raw(idx_ok), n, 'linear', 'extrap');
end

% ========================================================================
% 3. DETREND, FILTRADO Y RECORTE DE BORDES
% ========================================================================
% Usar un eje normalizado evita acondicionamiento innecesario en polyfit.
u = linspace(-1, 1, numel(phi_raw)).';
p = polyfit(u, phi_raw, 3);
phi_det = phi_raw - polyval(p, u);

% Filtrar SIEMPRE sobre la captura completa. Luego se recortan todos los
% canales con los mismos indices; asi no se crean nuevos transitorios al
% filtrar una senal que ya fue recortada.
x_hp_full = filtfilt(dHP, phi_det);
x_rr_full = filtfilt(dLP_rr, x_hp_full);

trim_s = max(1.5, min(trim_s, 0.12 * dur_cap_s));
n_trim = round(trim_s * Fs);
if numel(x_hp_full) > 2*n_trim + round(10*Fs)
    keep = (n_trim+1):(numel(x_hp_full)-n_trim);
else
    warning('Captura demasiado corta para recortar %.1f s por borde.', trim_s);
    keep = 1:numel(x_hp_full);
end

x = x_hp_full(keep);
x_rr = x_rr_full(keep);
phi_use = phi_raw(keep);

info.phi = phi_use(:);
info.x = x(:);
info.x_rr = x_rr(:);

fprintf('Senal util: %.1f s (%d muestras a %.1f Hz)\n', numel(x)/Fs, numel(x), Fs);

% ========================================================================
% 4. RESPIRACION
% ========================================================================
[f_rr_ax, P_rr, bin_hz] = espectro(x_rr, Fs, pad_mult);
info.f_rr_axis = f_rr_ax;
info.P_rr = P_rr;
mask_rr = (f_rr_ax >= rr_min_hz) & (f_rr_ax <= rr_max_hz);

[f_rr_est, rr_quality] = select_rr_peak_robust(f_rr_ax, P_rr, mask_rr, bin_hz);
info.rr_prominence = rr_quality.prominence;
info.rr_score_ratio = rr_quality.score_ratio;

if ~isfinite(f_rr_est) || f_rr_est < rr_min_hz || f_rr_est > rr_max_hz
    info.rr_fallback = true;
    rpm_est = NaN;
    info.f_rr_est = NaN;
else
    rpm_est = double(60*f_rr_est);
    info.f_rr_est = f_rr_est;
    info.rr_low_confidence = rr_quality.prominence < 2.0 || rr_quality.score_ratio < 1.10;
    if info.rr_low_confidence
        warning('RR %.1f rpm con baja confianza (prom %.2f, ratio %.2f).', ...
            rpm_est, rr_quality.prominence, rr_quality.score_ratio);
    end
end

% ========================================================================
% 5. FRECUENCIA CARDIACA
% ========================================================================
% La banda se fija independientemente de RR. Los armonicos de respiracion
% se tratan en el selector; no se desplaza artificialmente el limite bajo.
dHR = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', hr_min_hz, ...
    'HalfPowerFrequency2', hr_max_hz, ...
    'SampleRate', Fs);

x_hr_full = filtfilt(dHR, x_hp_full);
x_hr = x_hr_full(keep);
info.x_hr = x_hr(:);

[f_hr_ax, P_hr] = espectro_promediado(x_hr, Fs, welch_win_s, 0.25, pad_mult);
info.f_hr_axis = f_hr_ax;
info.P_hr = P_hr;

mask_hr = (f_hr_ax >= hr_min_hz) & (f_hr_ax <= hr_max_hz);
rr_for_hr = info.f_rr_est;

% Solo se elimina de forma dura la zona de los armonicos bajos. Los
% armonicos 4-6 se penalizan suavemente en el selector porque pueden caer
% exactamente sobre un HR real (ej.: 15 rpm -> 5o armonico = 75 bpm).
if isfinite(rr_for_hr)
    delta = max(0.05, 0.10*rr_for_hr);
    for h = 1:3
        mask_hr = mask_hr & ~(abs(f_hr_ax - h*rr_for_hr) <= delta);
    end
end

[f_hr_est, hr_quality] = select_hr_peak_robust(f_hr_ax, P_hr, mask_hr, rr_for_hr);
info.hr_prominence = hr_quality.prominence;
info.hr_score_ratio = hr_quality.score_ratio;

if ~isfinite(f_hr_est) || f_hr_est < hr_min_hz || f_hr_est > hr_max_hz
    info.hr_fallback = true;
    bpm_est = NaN;
    info.f_hr_est = NaN;
    info.hr_agreement = NaN;
    return;
end

bpm_est = double(60*f_hr_est);
info.f_hr_est = f_hr_est;

[info.hr_window_bpm, info.hr_agreement] = consistencia_hr( ...
    x_hr, Fs, hr_win_s, hr_min_hz, hr_max_hz, bpm_est, rr_for_hr);

n_hr_win = numel(info.hr_window_bpm);
info.hr_low_confidence = hr_quality.prominence < 2.0 || ...
                         hr_quality.score_ratio < 1.12 || ...
                         n_hr_win < 3 || ...
                         (isfinite(info.hr_agreement) && info.hr_agreement < 0.60);

fprintf('\n--- CONSISTENCIA HR ---\n');
fprintf('HR por ventana (bpm): %s\n', formatear_lista(info.hr_window_bpm));
if isfinite(info.hr_agreement)
    fprintf('Acuerdo con global (+-5 bpm): %.0f%%\n', 100*info.hr_agreement);
else
    fprintf('Acuerdo con global: no evaluable (muy pocas ventanas)\n');
end

if info.hr_low_confidence
    warning(['HR %.1f bpm con baja confianza (prom %.2f, ratio %.2f, ' ...
        '%d ventanas).'], bpm_est, hr_quality.prominence, ...
        hr_quality.score_ratio, n_hr_win);
end
end

% =========================================================================
% AJUSTE DE CIRCULO
% =========================================================================
function [c_t, c0, R0, n_fits] = estimar_centro_circulo(z, Fs, win_s)
z = z(:);
N = numel(z);
[c0, R0] = fit_circle_kasa(z);

Nw = round(win_s*Fs);
hop = max(1, round(Nw/4));
if Nw >= N || Nw < 32
    c_t = repmat(c0, N, 1);
    n_fits = 0;
    return;
end

centros = [];
tiempos = [];
for i0 = 1:hop:(N-Nw+1)
    seg = z(i0:i0+Nw-1);
    ang = unwrap(angle(seg-c0));
    if max(ang)-min(ang) < 0.35
        continue;
    end

    [cw, Rw] = fit_circle_kasa(seg);
    if ~isfinite(cw) || ~isfinite(Rw)
        continue;
    end
    if abs(cw-c0) > 0.6*R0 || Rw < 0.4*R0 || Rw > 2.5*R0
        continue;
    end

    centros(end+1,1) = cw; %#ok<AGROW>
    tiempos(end+1,1) = i0 + (Nw-1)/2; %#ok<AGROW>
end

n_fits = numel(centros);
if n_fits < 3
    c_t = repmat(c0, N, 1);
    return;
end

k_smooth = min(5, 2*floor(n_fits/4)+1);
centros = movmean(centros, k_smooth);
idx = (1:N).';
c_t = interp1(tiempos, centros, idx, 'linear');
c_t(idx < tiempos(1)) = centros(1);
c_t(idx > tiempos(end)) = centros(end);
end

function [c, R] = fit_circle_kasa(z)
x = real(z(:));
y = imag(z(:));
A = [x, y, ones(numel(x),1)];
b = x.^2 + y.^2;
p = A\b;
c = complex(p(1)/2, p(2)/2);
R = sqrt(max(p(3)+abs(c)^2, eps));
end

% =========================================================================
% CONSISTENCIA HR
% =========================================================================
function [bpm_win, acuerdo] = consistencia_hr(x_hr, Fs, win_s, fmin, fmax, bpm_ref, f_rr)
x_hr = x_hr(:);
N = numel(x_hr);
Nw = round(win_s*Fs);
hop = max(1, round(Nw/2));
bpm_win = [];

if Nw >= N
    acuerdo = NaN;
    return;
end

for i0 = 1:hop:(N-Nw+1)
    seg = x_hr(i0:i0+Nw-1);
    [f, P, ~] = espectro(seg, Fs, 8);
    m = (f >= fmin) & (f <= fmax);

    if isfinite(f_rr)
        delta = max(0.05, 0.10*f_rr);
        for h = 1:3
            m = m & ~(abs(f-h*f_rr) <= delta);
        end
    end

    [f0, ~] = select_hr_peak_robust(f, P, m, f_rr);
    if isfinite(f0)
        bpm_win(end+1,1) = 60*f0; %#ok<AGROW>
    end
end

if numel(bpm_win) < 2
    acuerdo = NaN;
else
    acuerdo = mean(abs(bpm_win-bpm_ref) <= 5);
end
end

% =========================================================================
% ESPECTROS
% =========================================================================
function [f_ax, P, bin_hz] = espectro(x, Fs, pad_mult, usar_hamming)
if nargin < 4, usar_hamming = true; end
x = x(:);
Nx = numel(x);
if Nx < 4
    f_ax = []; P = []; bin_hz = NaN; return;
end
Nfft = 2^nextpow2(max(4096, pad_mult*Nx));
bin_hz = Fs/Nfft;
if usar_hamming
    w = hamming(Nx);
else
    w = ones(Nx,1);
end
X = fft((x-mean(x)).*w, Nfft);
K = floor(Nfft/2)+1;
f_ax = (0:K-1).'*bin_hz;
P = abs(X(1:K));
end

function [f_ax, P] = espectro_promediado(x, Fs, win_s, hop_frac, pad_mult)
x = x(:);
N = numel(x);
L = round(win_s*Fs);
if L >= N
    [f_ax, P] = espectro(x, Fs, pad_mult);
    return;
end

hop = max(1, round(hop_frac*L));
w = hann(L);
w = w/sqrt(mean(w.^2));
inicios = 1:hop:(N-L+1);
Pmat = [];
for i = 1:numel(inicios)
    seg = x(inicios(i):inicios(i)+L-1);
    [f_ax, Pi, ~] = espectro(seg.*w, Fs, pad_mult, false);
    Pmat(:,end+1) = Pi(:); %#ok<AGROW>
end
P = mean(Pmat, 2);
end

% =========================================================================
% SELECTORES DE PICO
% =========================================================================
function [f_est, quality] = select_rr_peak_robust(f_ax, P, mask, bin_hz)
quality = struct('prominence', 0, 'score', 0, 'score_ratio', 0);
idx_band = find(mask);
if isempty(idx_band)
    f_est = NaN;
    return;
end

idx_peaks = local_peak_indices(P, idx_band);
if isempty(idx_peaks)
    [~, ir] = max(P(idx_band));
    idx_peaks = idx_band(ir);
end

scores = zeros(size(idx_peaks));
proms = zeros(size(idx_peaks));
f_cands = zeros(size(idx_peaks));
for i = 1:numel(idx_peaks)
    k = idx_peaks(i);
    f0 = parabolic_peak_frequency(f_ax, P, k, bin_hz);
    nf = local_noise_floor(f_ax, P, f0, 0.12, 0.025, mask);
    prom = P(k)/(nf+eps);

    edge_penalty = 1;
    if f0 < 8/60
        edge_penalty = 0.25;
    elseif f0 < 10/60
        edge_penalty = 0.55;
    elseif f0 < 12/60
        edge_penalty = 0.85;
    end

    scores(i) = prom*edge_penalty;
    proms(i) = prom;
    f_cands(i) = f0;
end

scores(~isfinite(scores)) = 0;
if ~any(scores > 0)
    f_est = NaN;
    return;
end

[~, best] = max(scores);
f_est = f_cands(best);
f_est = prefer_rr_fundamental(f_est, f_ax, P, mask, bin_hz);

k_final = nearest_masked_bin(f_ax, mask, f_est);
if ~isempty(k_final)
    nf = local_noise_floor(f_ax, P, f_est, 0.12, 0.025, mask);
    quality.prominence = P(k_final)/(nf+eps);
end
quality.score = scores(best);
sorted_scores = sort(scores(:), 'descend');
if numel(sorted_scores) >= 2 && sorted_scores(2) > 0
    quality.score_ratio = sorted_scores(1)/sorted_scores(2);
else
    quality.score_ratio = Inf;
end
end

function f_est = prefer_rr_fundamental(f_est, f_ax, P, mask, bin_hz)
if ~isfinite(f_est) || f_est < 24/60
    return;
end

k_est = nearest_masked_bin(f_ax, mask, f_est);
if isempty(k_est), return; end
p_est = P(k_est);
best_f = f_est;
best_ratio = 0;

for divisor = [2 3]
    f_candidate = f_est/divisor;
    if f_candidate < 8/60 || f_candidate > 24/60
        continue;
    end
    if divisor == 3 && (f_est < 33/60 || f_candidate < 10/60)
        continue;
    end

    k0 = nearest_masked_bin(f_ax, mask, f_candidate);
    if isempty(k0), continue; end
    half_width_bins = max(2, round(0.035/bin_hz));
    lo = max(1, k0-half_width_bins);
    hi = min(numel(P), k0+half_width_bins);
    local_idx = lo:hi;
    local_idx = local_idx(mask(local_idx));
    if isempty(local_idx), continue; end

    [p_fund, ir] = max(P(local_idx));
    f_fund = f_ax(local_idx(ir));
    ratio = p_fund/(p_est+eps);
    required = 0.45;
    if divisor == 3, required = 0.55; end

    if ratio >= required && ratio > best_ratio
        best_ratio = ratio;
        best_f = f_fund;
    end
end

f_est = best_f;
end

function [f_est, quality] = select_hr_peak_robust(f_ax, P, mask, f_rr)
quality = struct('prominence', 0, 'score', 0, 'score_ratio', 0);
idx_band = find(mask);
if isempty(idx_band)
    f_est = NaN;
    return;
end

idx_peaks = local_peak_indices(P, idx_band);
if isempty(idx_peaks)
    [~, ir] = max(P(idx_band));
    idx_peaks = idx_band(ir);
end

scores = zeros(size(idx_peaks));
proms = zeros(size(idx_peaks));
f_cands = zeros(size(idx_peaks));

for i = 1:numel(idx_peaks)
    k = idx_peaks(i);
    f0 = f_ax(k);
    nf = local_noise_floor(f_ax, P, f0, 0.35, 0.06, mask);
    prom = P(k)/(nf+eps);
    proms(i) = prom;

    harm_pen = 0;
    if isfinite(f_rr) && f_rr > 0
        delta_low = max(0.05, 0.10*f_rr);
        delta_high = max(0.035, 0.07*f_rr);
        for h = 1:3
            harm_pen = harm_pen + exp(-((f0-h*f_rr).^2)/(2*delta_low^2));
        end
        for h = 4:6
            harm_pen = harm_pen + 0.12*exp(-((f0-h*f_rr).^2)/(2*delta_high^2));
        end
    end

    adult_bonus = 1;
    if f0 >= 60/60 && f0 <= 130/60
        adult_bonus = 1.15;
    end

    high_rate_penalty = 1;
    if f0 > 130/60, high_rate_penalty = 0.70; end
    if f0 > 150/60, high_rate_penalty = 0.55; end

    scores(i) = prom*adult_bonus*high_rate_penalty/(1+harm_pen);
    f_cands(i) = f0;
end

scores(~isfinite(scores)) = 0;
if ~any(scores > 0)
    f_est = NaN;
    return;
end

[~, best] = max(scores);
f_est = f_cands(best);
quality.prominence = proms(best);
quality.score = scores(best);
sorted_scores = sort(scores(:), 'descend');
if numel(sorted_scores) >= 2 && sorted_scores(2) > 0
    quality.score_ratio = sorted_scores(1)/sorted_scores(2);
else
    quality.score_ratio = Inf;
end
end

function idx_peaks = local_peak_indices(P, idx_band)
idx_peaks = [];
for ii = 1:numel(idx_band)
    k = idx_band(ii);
    if k > 1 && k < numel(P) && P(k) >= P(k-1) && P(k) >= P(k+1)
        idx_peaks(end+1) = k; %#ok<AGROW>
    end
end
end

function f0 = parabolic_peak_frequency(f_ax, P, k, bin_hz)
f0 = f_ax(k);
if k <= 1 || k >= numel(P), return; end
al = P(k-1); be = P(k); ga = P(k+1);
denom = al - 2*be + ga;
if abs(denom) > 1e-12
    di = 0.5*(al-ga)/denom;
    di = max(-0.5, min(0.5, di));
    f0 = f_ax(k) + di*bin_hz;
end
end

function nf = local_noise_floor(f_ax, P, f0, half_width, guard_width, mask)
idx = mask & abs(f_ax-f0) <= half_width & abs(f_ax-f0) >= guard_width;
if nnz(idx) < 8
    idx = mask;
end
if ~any(idx)
    nf = NaN;
    return;
end
nf = median(P(idx));
if ~isfinite(nf) || nf <= 0
    nf = median(P(mask)) + eps;
end
end

function k = nearest_masked_bin(f_ax, mask, f0)
idx = find(mask);
if isempty(idx) || ~isfinite(f0)
    k = [];
    return;
end
[~, ir] = min(abs(f_ax(idx)-f0));
k = idx(ir);
end

function txt = formatear_lista(v)
if isempty(v)
    txt = 'sin datos';
else
    parts = strings(1, numel(v));
    for i = 1:numel(v)
        parts(i) = sprintf('%.1f', v(i));
    end
    txt = strjoin(parts, ', ');
end
end

function check_filter_fs(d, name, fs_runtime)
try
    fs_design = d.SampleRate;
catch
    return;
end
if abs(fs_design-fs_runtime) > 1e-6
    error('%s fue disenado para %.3f Hz, pero Fs actual es %.3f Hz.', ...
        name, fs_design, fs_runtime);
end
end

% =========================================================================
% OFFSET RESIDUAL DE PORTADORA
% =========================================================================
function [z_out, f_off, vueltas, deriva_pp] = quitar_offset_portadora(z, Fs, fc_track)
if nargin < 3 || isempty(fc_track), fc_track = 0.06; end

N = numel(z);
n = (0:N-1).';

% Pico grueso cerca de DC.
Nf = 2^nextpow2(8*N);
Zf = fft(z, Nf);
f = (0:Nf-1).'/Nf*Fs;
f(f > Fs/2) = f(f > Fs/2)-Fs;
mag = abs(Zf);
mag(abs(f) > 2) = 0;
[~, k] = max(mag);
f0 = f(k);

z1 = z .* exp(-1j*2*pi*f0*n/Fs);

% Refinamiento mediante pendiente lineal de fase.
psi = unwrap(angle(z1));
u = linspace(-1,1,N).';
p = polyfit(u, psi, 1);
% du/dn = 2/(N-1). Convertir pendiente por u a rad/muestra.
rad_per_sample = p(1) * 2/max(N-1,1);
f1 = rad_per_sample*Fs/(2*pi);

z2 = z1 .* exp(-1j*2*pi*f1*n/Fs);
f_off = f0 + f1;
vueltas = f_off*N/Fs;

if fc_track <= 0
    z_out = z2;
    deriva_pp = 0;
    return;
end

psi = unwrap(angle(z2));
Wn = fc_track/(Fs/2);
[bb, aa] = butter(2, Wn, 'low');
psi_lento = filtfilt(bb, aa, psi);
z_out = z2 .* exp(-1j*psi_lento);
deriva_pp = max(psi_lento)-min(psi_lento);
end
