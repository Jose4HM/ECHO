function [rpm_est, bpm_est, info] = procesar_signos_vitales_fixed(buffer_IQ, Fs, dHP, dLP_rr)
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
%   info    – Struct con banderas de diagnostico (para no confundir un
%             valor generico con una deteccion real):
%               .rr_fallback       true si RR no encontro pico real en
%                                   banda y se devolvio el generico (15 rpm)
%               .hr_fallback       true si HR no encontro ningun
%                                   candidato en banda y se devolvio el
%                                   generico (75 bpm)
%               .hr_low_confidence true si HR si encontro un pico pero
%                                   con poca prominencia/ratio de confianza
%               .hr_prominence, .hr_score_ratio  metricas crudas de HR

pad_mult = 8;   % Zero-padding para interpolación espectral fina

% Banderas de diagnostico: por defecto asumimos deteccion real; se
% marcan como fallback solo en los puntos donde el codigo realmente
% recurre a un valor fisiologico generico en vez de una medicion.
info = struct('rr_fallback', false, 'hr_fallback', false, ...
    'hr_low_confidence', false, 'hr_prominence', 0, 'hr_score_ratio', 0);

% Validar que los filtros pre-diseñados (dHP, dLP_rr) sean consistentes
% con el Fs de este buffer. Si el objeto de filtro expone 'SampleRate'
% (p. ej. digitalFilter de designfilt), se compara; si no lo expone, se
% omite la verificación silenciosamente (no todos los tipos de filtro
% de MATLAB exponen esta propiedad de forma uniforme). Un desajuste aquí
% no lanza error en filtfilt -- simplemente filtra en banda equivocada
% sin avisar, así que vale la pena advertir cuando se puede detectar.
if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('Fs debe ser un escalar positivo en Hz.');
end

buffer_IQ = buffer_IQ(:);

check_filter_fs = @(d, name) local_check_filter_fs(d, name, Fs);
check_filter_fs(dHP, 'dHP');
check_filter_fs(dLP_rr, 'dLP_rr');

% ================================================================
% 1. ACONDICIONAMIENTO IQ
% ================================================================
% Validar entrada: NaN/Inf contaminan silenciosamente toda la cadena
% (mean, std, filtfilt) sin lanzar ningún error visible.
if any(~isfinite(buffer_IQ))
    n_bad = sum(~isfinite(buffer_IQ));
    warning('buffer_IQ contiene %d muestras no finitas (NaN/Inf); se interpolarán.', n_bad);
    idx_good = isfinite(buffer_IQ);
    if sum(idx_good) < 2
        error('buffer_IQ no tiene suficientes muestras finitas para procesar.');
    end
    idx_all = (1:numel(buffer_IQ)).';
    buffer_IQ = interp1(idx_all(idx_good), buffer_IQ(idx_good), idx_all, 'linear', 'extrap');
end

N = length(buffer_IQ);

if N < 4 * Fs
    warning('buffer_IQ tiene %.2f s de datos; BPM/RPM seran poco estables.', N / Fs);
end

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

% La fase es poco confiable cuando la trayectoria IQ pasa muy cerca del
% origen. En capturas reales eso mete saltos que despues pueden aparecer
% como picos falsos en la banda cardiaca. Se reemplazan solo esas muestras
% por interpolacion temporal de la fase ya desenvuelta.
iq_radius = abs(S_iq);
radius_floor = max(prctile(iq_radius, 10), 0.20 * median(iq_radius));
idx_phase_good = isfinite(phi_raw) & iq_radius > radius_floor;
if nnz(idx_phase_good) >= max(16, round(0.25 * numel(phi_raw)))
    phi_raw = interp1(n(idx_phase_good), phi_raw(idx_phase_good), n, ...
        'linear', 'extrap');
else
    warning('Muchas muestras IQ tienen baja magnitud; la fase puede ser inestable.');
end

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

% Banda fisiologica adulta de respiracion: 8-35 RPM
rr_min_hz = 6 / 60;
rr_max_hz = 45 / 60;
mask_rr  = (f_rr_ax >= rr_min_hz) & (f_rr_ax <= rr_max_hz);
f_rr_est = NaN;
legacy_rr_selector_enabled = ~isempty(getenv('VITAL_USE_LEGACY_RR'));

if legacy_rr_selector_enabled && any(mask_rr)
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

% En las pruebas reales esto sesga 24-26 RPM hacia 6-10 RPM.
f_rr_est = f_peak;
end

% Fallback si no se encontró nada razonable
if any(mask_rr)
    f_rr_est = select_rr_peak_robust(f_rr_ax, P_rr, mask_rr, bin_hz);
end

if isnan(f_rr_est) || f_rr_est < rr_min_hz || f_rr_est > rr_max_hz
    info.rr_fallback = true;   % No se encontro un pico real en banda
    f_rr_est = 0.25;   % 15 RPM: valor adulto en reposo típico
end
rpm_est = double(60 * f_rr_est(1));   % Garantizar escalar

% ================================================================
% 4. CANCELACIÓN DE RESPIRACIÓN (SIC — Successive Interference Cancellation)
% ================================================================
% ✅ FIX BUG 1: construir modelo sobre x_rr Y sustraer de x_rr → residuo
% coherente. Luego restar de x (sin filtrar) solo si x_rr ≈ x en banda HR.
tbb = (0:numel(x)-1).' / Fs;
% Normalizar la rampa temporal (centrada y escalada a [-1,1] aprox.) para
% evitar mal condicionamiento de H: sin esto, la columna de rampa lineal
% puede tener magnitud mucho mayor que las columnas senoidales, lo que
% degrada la estabilidad numérica de la regresión por mínimos cuadrados.
tbb_n = (tbb - mean(tbb)) / (max(tbb) - min(tbb) + eps);
H   = [ones(numel(tbb),1), tbb_n, ...
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
fmin_hr = max(45/60, f_rr_est + margin);
fmax_hr = 3.33;   % 200 BPM máximo razonable
if fmin_hr >= fmax_hr - 0.1
    fmin_hr = 45/60;
end




% Bandpass cardiaco estable con SampleRate explicito.
dHR = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', fmin_hr, ...
    'HalfPowerFrequency2', fmax_hr, ...
    'SampleRate', Fs);
x_hr = filtfilt(dHR, x_rrfree);

Nfft_hr  = 2^nextpow2(pad_mult * numel(x_hr));
bin_hr   = Fs / Nfft_hr;
X_hr     = fft(x_hr .* hamming(numel(x_hr)), Nfft_hr);
K_hr     = floor(Nfft_hr/2) + 1;
f_hr_ax  = (0:K_hr-1) * bin_hr;
P_hr     = abs(X_hr(1:K_hr));


% % --- INICIO BLOQUE DE DEBUG ---
% % x_hr es la señal cardíaca filtrada y lista para la FFT
% N_debug = numel(x_hr);
% f_debug = (0:N_debug-1) * (Fs / N_debug);
% X_debug = abs(fft(x_hr .* hamming(N_debug), N_debug)); 
% 
% subplot(2, 2, [3, 4]);
% % figure;
% plot(f_debug(1:floor(N_debug/2)), X_debug(1:floor(N_debug/2)));
% title('Espectro de frecuencia de la señal cardíaca (x_hr)');
% xlabel('Frecuencia (Hz)');
% ylabel('Magnitud');
% grid on;
% % --- FIN BLOQUE DE DEBUG ---

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
    info.hr_fallback = true;   % No hay ningun candidato en banda
    bpm_est = 75;   % Fallback
    return;
end

% Penalización gaussiana sobre armónicos de respiración
pen = zeros(size(cands_f));
for h = 1:6
    pen = pen + exp(-((cands_f - h*f_rr_est).^2) / (2 * delta^2));
end
score = cands_P ./ (1 + 3 * pen);

% Guardia contra NaN/Inf en score (división por cero o cands_P corruptos)
score(~isfinite(score)) = 0;

if ~any(score > 0)
    % Sin candidato válido: devolver fallback fisiológico
    info.hr_fallback = true;
    bpm_est = 75;
    return;
end

% max() garantiza escalar; ii(1) cubre el caso raro de empate exacto
[~, ii] = max(score);
ii      = ii(1);
bpm_est = double(60 * cands_f(ii));
hr_quality = struct('prominence', 0, 'score_ratio', 0);
if any(mask_final)
    [f_hr_est, hr_quality] = select_hr_peak_robust(f_hr_ax, P_hr, ...
        mask_final, f_rr_est);
    bpm_est = 60 * f_hr_est;
end

info.hr_prominence = hr_quality.prominence;
info.hr_score_ratio = hr_quality.score_ratio;
info.hr_low_confidence = hr_quality.prominence < 2.0 || hr_quality.score_ratio < 1.12;

% Guardia fisiológica final — también garantiza escalar
if info.hr_low_confidence
    warning(['HR con baja confianza: prominencia %.2f, ratio %.2f. ' ...
        'Si no coincide con el valor real, probablemente domina residuo ' ...
        'respiratorio/multipath.'], ...
        hr_quality.prominence, hr_quality.score_ratio);
end

bpm_est = double(max(30, min(200, bpm_est(1))));
end

% ================================================================
% FUNCIÓN AUXILIAR LOCAL
% ================================================================
function f_est = select_rr_peak_robust(f_ax, P, mask, bin_hz)
idx_band = find(mask);
if isempty(idx_band)
    f_est = NaN;
    return;
end

idx_peaks = local_peak_indices(P, idx_band);
if isempty(idx_peaks)
    [~, i_rel] = max(P(idx_band));
    idx_peaks = idx_band(i_rel);
end

scores = zeros(size(idx_peaks));
f_cands = zeros(size(idx_peaks));
for i = 1:numel(idx_peaks)
    k = idx_peaks(i);
    f0 = parabolic_peak_frequency(f_ax, P, k, bin_hz);
    local_floor = local_noise_floor(f_ax, P, f0, 0.12, 0.025, mask);
    prom = P(k) / (local_floor + eps);

    % Penalizacion del borde inferior: en capturas reales, el drift lento
    % puede competir con respiraciones de 12-18 rpm y arrastrar el resultado
    % hacia 6-9 rpm. No se bloquea, pero debe ser muy prominente para ganar.
    edge_penalty = 1.0;
    if f0 < 8/60
        edge_penalty = 0.25;
    elseif f0 < 10/60
        edge_penalty = 0.55;
    elseif f0 < 12/60
        edge_penalty = 0.85;
    end

    scores(i) = prom * edge_penalty;
    f_cands(i) = f0;
end

scores(~isfinite(scores)) = 0;
if ~any(scores > 0)
    [~, i_rel] = max(P(idx_band));
    f_est = f_ax(idx_band(i_rel));
else
    [~, best] = max(scores);
    f_est = f_cands(best);
end

f_est = prefer_rr_fundamental(f_est, f_ax, P, mask, bin_hz);
end

function f_est = prefer_rr_fundamental(f_est, f_ax, P, mask, bin_hz)
% Si el pico de respiracion elegido parece ser un armonico, preferir la
% fundamental cuando conserva energia suficiente. Esto corrige casos como
% 28-30 rpm cuando la respiracion real esta cerca de 14-15 rpm.
if ~isfinite(f_est) || f_est < 24/60
    return;
end

k_est = nearest_masked_bin(f_ax, mask, f_est);
if isempty(k_est)
    return;
end
p_est = P(k_est);

best_f = f_est;
best_ratio = 0;
for divisor = [2 3]
    f_candidate = f_est / divisor;
    if f_candidate < 8/60 || f_candidate > 24/60
        continue;
    end

    % Dividir por 3 es mas riesgoso: 27 rpm / 3 = 9 rpm suele ser drift,
    % no respiracion real. Solo se permite cuando el pico alto es compatible
    % con un tercer armonico claro y la fundamental resultante no cae al borde.
    if divisor == 3 && (f_est < 33/60 || f_candidate < 10/60)
        continue;
    end

    k0 = nearest_masked_bin(f_ax, mask, f_candidate);
    if isempty(k0)
        continue;
    end

    half_width_bins = max(2, round(0.035 / bin_hz));
    lo = max(1, k0 - half_width_bins);
    hi = min(numel(P), k0 + half_width_bins);
    local_idx = lo:hi;
    local_idx = local_idx(mask(local_idx));
    if isempty(local_idx)
        continue;
    end

    [p_fund, i_rel] = max(P(local_idx));
    f_fund = f_ax(local_idx(i_rel));
    ratio = p_fund / (p_est + eps);
    required_ratio = 0.45;
    if divisor == 3
        required_ratio = 0.55;
    end

    if ratio >= required_ratio && ratio > best_ratio
        best_ratio = ratio;
        best_f = f_fund;
    end
end

if best_ratio > 0
    f_est = best_f;
end
end

function k = nearest_masked_bin(f_ax, mask, f0)
idx = find(mask);
if isempty(idx)
    k = [];
    return;
end
[~, i_rel] = min(abs(f_ax(idx) - f0));
k = idx(i_rel);
end

function [f_est, quality] = select_hr_peak_robust(f_ax, P, mask, f_rr)
quality = struct('prominence', 0, 'score_ratio', 0);
idx_band = find(mask);
if isempty(idx_band)
    f_est = 75 / 60;
    return;
end

idx_peaks = local_peak_indices(P, idx_band);
if isempty(idx_peaks)
    [~, i_rel] = max(P(idx_band));
    idx_peaks = idx_band(i_rel);
end

scores = zeros(size(idx_peaks));
f_cands = zeros(size(idx_peaks));
for i = 1:numel(idx_peaks)
    k = idx_peaks(i);
    f0 = f_ax(k);
    local_floor = local_noise_floor(f_ax, P, f0, 0.35, 0.06, mask);
    prom = P(k) / (local_floor + eps);

    harm_pen = 0;
    delta_low = max(0.05, 0.10 * f_rr);
    delta_high = max(0.035, 0.07 * f_rr);
    for h = 1:3
        harm_pen = harm_pen + 1.0 * ...
            exp(-((f0 - h*f_rr).^2) / (2 * delta_low^2));
    end
    for h = 4:6
        harm_pen = harm_pen + 0.12 * ...
            exp(-((f0 - h*f_rr).^2) / (2 * delta_high^2));
    end

    % Priorizacion suave del rango adulto comun, sin excluir bradicardia
    % ni taquicardia si tienen mejor prominencia.
    adult_bonus = 1.0;
    if f0 >= 60/60 && f0 <= 130/60
        adult_bonus = 1.15;
    end

    high_rate_penalty = 1.0;
    if f0 > 130/60
        high_rate_penalty = 0.70;
    end
    if f0 > 150/60
        high_rate_penalty = 0.55;
    end

    scores(i) = prom * adult_bonus * high_rate_penalty / (1 + 1 * harm_pen);
    f_cands(i) = f0;
end

scores(~isfinite(scores)) = 0;
if ~any(scores > 0)
    [~, i_rel] = max(P(idx_band));
    f_est = f_ax(idx_band(i_rel));
else
    [~, best] = max(scores);
    f_est = f_cands(best);
    quality.prominence = scores(best);
    sorted_scores = sort(scores(:), 'descend');
    if numel(sorted_scores) >= 2 && sorted_scores(2) > 0
        quality.score_ratio = sorted_scores(1) / sorted_scores(2);
    else
        quality.score_ratio = Inf;
    end
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
if k <= 1 || k >= numel(P)
    return;
end
al = P(k-1);
be = P(k);
ga = P(k+1);
denom = al - 2*be + ga;
if abs(denom) > 1e-12
    delta_idx = 0.5 * (al - ga) / denom;
    delta_idx = max(-0.5, min(0.5, delta_idx));
    f0 = f_ax(k) + delta_idx * bin_hz;
end
end

function nf = local_noise_floor(f_ax, P, f0, half_width, guard_width, mask)
idx = mask & abs(f_ax - f0) <= half_width & abs(f_ax - f0) >= guard_width;
if nnz(idx) < 8
    idx = mask;
end
nf = median(P(idx));
if ~isfinite(nf) || nf <= 0
    nf = median(P(mask)) + eps;
end
end

function local_check_filter_fs(d, name, fs_runtime)
% Compara el Fs de diseño del filtro (si está expuesto) contra el Fs
% real con el que se está llamando esta función. Un desajuste indica
% que el filtro fue diseñado para otra tasa de muestreo y va a filtrar
% en la banda de frecuencia equivocada sin lanzar ningún error.
try
    fs_design = d.SampleRate; 
catch
    return;  % El objeto no expone SampleRate; no se puede verificar.
end
if abs(fs_design - fs_runtime) > 1e-6
    warning(['%s fue diseñado para Fs=%.3f Hz pero se está llamando con ' ...
        'Fs=%.3f Hz. El filtro actuará sobre la banda de frecuencia ' ...
        'equivocada sin lanzar error -- revisa el diseño del filtro.'], ...
        name, fs_design, fs_runtime);
end
end
