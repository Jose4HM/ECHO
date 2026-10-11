% =========================================================================
% MASTER SCRIPT: RADAR DE SIGNOS VITALES
% =========================================================================
clear; clc; close all;

% --- PASO 1: CAPTURA ---
try
    [rx_continuous, cfg] = cw();
catch ME
    error(['Fallo en cw.m durante la captura: %s\n' ...
        'Revisa conexion del hardware y configuracion antes de reintentar.'], ...
        ME.message);
end

if ~exist('rx_continuous', 'var') || isempty(rx_continuous)
    error('cw.m no genero rx_continuous.');
end

fprintf('\nEjecutando Front-End DSP...\n');
[iq_frontend, Fs_frontend] = frontend_dsp(rx_continuous, cfg);

if isempty(iq_frontend) || ~isfinite(Fs_frontend) || Fs_frontend <= 0
    error('frontend_dsp devolvio una salida invalida.');
end

% --- PASO 2: FILTROS FISIOLOGICOS ---
dHP = designfilt('highpassiir', 'FilterOrder', 2, ...
    'HalfPowerFrequency', 0.1, 'SampleRate', Fs_frontend);
 

dLP_rr = designfilt('lowpassiir', 'FilterOrder', 2, ...
    'HalfPowerFrequency', 0.8, 'SampleRate', Fs_frontend);

% --- PASO 3: PROCESAMIENTO ---
fprintf('Datos capturados. Procesando signos vitales...\n');
[rpm, bpm, vital_info] = procesar_signos_vitales_fixed(iq_frontend, Fs_frontend, dHP, dLP_rr);

% --- PASO 4: RESULTADOS ---
fprintf('=============================\n');
fprintf('RESULTADOS FINALES:\n');
fprintf('Respiracion: %.1f rpm\n', rpm);
fprintf('Ritmo cardiaco: %.1f bpm\n', bpm);
fprintf('=============================\n');

if vital_info.rr_fallback
    fprintf(2, ['AVISO: RR = %.1f rpm es un valor GENERICO de fallback, ' ...
        'no una deteccion real. No hubo ningun pico valido en la banda ' ...
        'de respiracion.\n'], rpm);
end
if vital_info.hr_fallback
    fprintf(2, ['AVISO: HR = %.1f bpm es un valor GENERICO de fallback, ' ...
        'no una deteccion real. No hubo ningun candidato en la banda ' ...
        'cardiaca.\n'], bpm);
elseif vital_info.hr_low_confidence
    fprintf(2, ['AVISO: HR = %.1f bpm proviene de un pico real pero de ' ...
        'baja confianza (prominencia %.2f, ratio %.2f). Tratar con ' ...
        'cautela.\n'], bpm, vital_info.hr_prominence, vital_info.hr_score_ratio);
end

% --- PASO 5: GRAFICAS DE DIAGNOSTICO ---
plot_diagnostico_radar(rx_continuous, cfg, iq_frontend, Fs_frontend, ...
    dHP, dLP_rr, rpm, bpm);

function plot_diagnostico_radar(rx_raw, cfg, iq_bb, Fs_bb, dHP, dLP_rr, rpm, bpm)
%PLOT_DIAGNOSTICO_RADAR Resume el flujo de procesamiento en tiempo/frecuencia.

Fs_raw = cfg.fs;
rx_raw = rx_raw(:);
iq_bb = iq_bb(:);

diag = calcular_diagnostico_vital(iq_bb, Fs_bb, dHP, dLP_rr, rpm, bpm);
imprimir_candidatos_diagnostico(diag);

% Creamos la figura SIN forzar el color de fondo (deja que MATLAB use el tema actual)
figure('Name', 'Diagnostico radar signos vitales', ...
    'Units', 'normalized', ...
    'Position', [0.05 0.08 0.90 0.82]);

% El título general
sgtitle(sprintf('Diagnostico radar | RR %.1f rpm | HR %.1f bpm', rpm, bpm), ...
    'FontWeight', 'bold', 'FontSize', 14); % Ajuste de tamaño para que se vea bien

% 1) Envolvente de amplitud RX por bloque, a lo largo de TODA la
% captura (no solo los primeros ms). Un punto por cada intento de
% transmit+receive, tomado directamente de cw.m, asi que no depende del
% indexado interno de rx_continuous ni se ve afectado si algun bloque
% salio incompleto. Los bloques incompletos se marcan aparte para poder
% distinguir "hubo señal real y baja" de "el bloque en si fallo".
subplot(3, 2, 1);
if isfield(cfg, 'block_peak') && isfield(cfg, 'block_rms')
    n_blk = numel(cfg.block_peak);
    t_blk_s = (0:n_blk-1) * (cfg.block_len / cfg.fs);
    plot(t_blk_s, cfg.block_peak, 'b'); hold on;
    plot(t_blk_s, cfg.block_rms, 'c');
    yline(0.99, ':r');
    idx_incompl = cfg.actual_counts ~= cfg.block_len;
    leg = {'pico |bloque|', 'RMS bloque', 'umbral saturacion'};
    if any(idx_incompl)
        plot(t_blk_s(idx_incompl), cfg.block_peak(idx_incompl), 'rx', ...
            'MarkerSize', 6, 'LineWidth', 1.2);
        leg{end+1} = 'bloque incompleto';
    end
    grid on;
    xlabel('s');
    ylabel('amplitud');
    title(sprintf('Envolvente RX, captura completa (%d/%d bloques incompletos)', ...
        nnz(idx_incompl), n_blk));
    legend(leg, 'Location', 'best');
else
    % cw.m antiguo sin instrumentacion por bloque: mostrar solo el
    % primer tramo, como antes.
    n_raw = min(numel(rx_raw), round(0.01 * Fs_raw));
    t_raw_ms = (0:n_raw-1) / Fs_raw * 1e3;
    plot(t_raw_ms, real(rx_raw(1:n_raw)), 'b');
    hold on;
    plot(t_raw_ms, imag(rx_raw(1:n_raw)), 'r');
    grid on;
    xlabel('ms');
    ylabel('amplitud');
    title('Senal recibida RX @1 MHz (solo primeros ms; cw.m sin instrumentar)');
    legend('I', 'Q', 'Location', 'best');
end

% 2) Constelacion IQ baseband.
subplot(3, 2, 2);
idx_const = round(linspace(1, numel(iq_bb), min(6000, numel(iq_bb))));
plot(real(iq_bb(idx_const)), imag(iq_bb(idx_const)), '.', 'MarkerSize', 3);
axis equal;
grid on;
xlabel('I');
ylabel('Q');
title(sprintf('Constelacion IQ frontend @%.0f Hz', Fs_bb));

% 3) FFT de la senal cruda: debe verse el tono low-IF cerca de 100 kHz.
subplot(3, 2, 3);
plot_fft_complex(rx_raw, Fs_raw, 262144, [-150 150], 1e3);
xline(100, '--r', '100 kHz', 'LabelVerticalAlignment', 'bottom');
grid on;
xlabel('kHz');
ylabel('dB');
title('FFT RX cruda antes del frontend');

% 4) FFT tras frontend: el contenido util ya esta cerca de 0 Hz.
subplot(3, 2, 4);
plot_fft_complex(iq_bb, Fs_bb, 65536, [-150 150], 1);
grid on;
xlabel('Hz');
ylabel('dB');
title('FFT IQ baseband despues del frontend');

% 5) FFT de respiracion despues de HP + LP.
subplot(3, 2, 5);
plot(diag.f_rr, diag.P_rr_db, 'LineWidth', 1.2);
hold on;
plot_candidatos_fft(diag.rr_cand_f, diag.rr_cand_db, 60, 'rpm');
xline(diag.f_rr_est, '--r', sprintf('%.1f rpm', rpm), ...
    'LabelVerticalAlignment', 'bottom');
xline(0.10, ':k', '6 rpm');
xline(0.75, ':k', '45 rpm');
grid on;
xlim([0 1.2]);
xlabel('Hz');
ylabel('dB');
title('FFT respiracion: fase -> HP -> LP');

% 6) FFT cardiaca despues de SIC + bandpass HR.
subplot(3, 2, 6);
plot(diag.f_hr, diag.P_hr_db, 'LineWidth', 1.2);
hold on; 
plot_candidatos_fft(diag.hr_cand_f, diag.hr_cand_db, 60, 'bpm');
xline(diag.f_hr_est, '--r', sprintf('%.1f bpm', bpm), ...
    'LabelVerticalAlignment', 'bottom');
xline(45/60, ':k', '45 bpm');
xline(3.33, ':k', '200 bpm');
grid on;
xlim([0.4 3.6]);
xlabel('Hz');
ylabel('dB');
title('FFT corazon: SIC -> bandpass HR');
end

function diag = calcular_diagnostico_vital(iq_bb, Fs, dHP, dLP_rr, rpm, bpm)
I = real(iq_bb) - mean(real(iq_bb));
Q = imag(iq_bb) - mean(imag(iq_bb));

alpha = dot(I, Q) / (dot(I, I) + eps);
Q_ortho = Q - alpha * I;
g = std(I) / (std(Q_ortho) + eps);
S_iq = I + 1j * (Q_ortho * g);

phi_raw = unwrap(angle(S_iq));
n = (0:numel(phi_raw)-1).';

iq_radius = abs(S_iq);
radius_floor = max(prctile(iq_radius, 10), 0.20 * median(iq_radius));
idx_phase_good = isfinite(phi_raw) & iq_radius > radius_floor;
if nnz(idx_phase_good) >= max(16, round(0.25 * numel(phi_raw)))
    phi_raw = interp1(n(idx_phase_good), phi_raw(idx_phase_good), n, ...
        'linear', 'extrap');
end

p = polyfit(n, phi_raw, 1);
phi_cfo = phi_raw - polyval(p, n);
x = filtfilt(dHP, phi_cfo);
x_rr = filtfilt(dLP_rr, x);

f_rr_est = max(6/60, min(45/60, rpm / 60));
t = (0:numel(x)-1).' / Fs;
t_n = (t - mean(t)) / (max(t) - min(t) + eps);
H = [ones(numel(t), 1), t_n, ...
    sin(2*pi*  f_rr_est*t), cos(2*pi*  f_rr_est*t), ...
    sin(2*pi*2*f_rr_est*t), cos(2*pi*2*f_rr_est*t), ...
    sin(2*pi*3*f_rr_est*t), cos(2*pi*3*f_rr_est*t)];
theta = H \ x_rr;
x_rrfree = x - H * theta;

margin = max(0.04, 0.08 * f_rr_est);
fmin_hr = max(45/60, f_rr_est + margin);
fmax_hr = 3.33;
if fmin_hr >= fmax_hr - 0.1
    fmin_hr = 45/60;
end
dHR = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', fmin_hr, ...
    'HalfPowerFrequency2', fmax_hr, ...
    'SampleRate', Fs);
x_hr = filtfilt(dHR, x_rrfree);

[diag.f_rr, diag.P_rr_db] = single_sided_fft_db(x_rr, Fs, 8);
[diag.f_hr, diag.P_hr_db] = single_sided_fft_db(x_hr, Fs, 8);
diag.f_rr_est = f_rr_est;
diag.f_hr_est = max(0.5, min(3.5, bpm / 60));
[diag.rr_cand_f, diag.rr_cand_db] = top_spectral_candidates( ...
    diag.f_rr, diag.P_rr_db, diag.f_rr >= 6/60 & diag.f_rr <= 45/60, ...
    5, 0.035);
[diag.hr_cand_f, diag.hr_cand_db] = top_spectral_candidates( ...
    diag.f_hr, diag.P_hr_db, diag.f_hr >= 45/60 & diag.f_hr <= 3.33, ...
    5, 0.08);
end

function imprimir_candidatos_diagnostico(diag)
fprintf('\n--- CANDIDATOS DIAGNOSTICO ---\n');
fprintf('RR candidatos (rpm): %s\n', formatear_lista(diag.rr_cand_f * 60));
fprintf('HR candidatos (bpm): %s\n', formatear_lista(diag.hr_cand_f * 60));
fprintf('------------------------------\n');
end

function txt = formatear_lista(x)
if isempty(x)
    txt = 'sin candidatos';
    return;
end
parts = strings(1, numel(x));
for i = 1:numel(x)
    parts(i) = sprintf('%.1f', x(i));
end
txt = strjoin(parts, ', ');
end

function plot_candidatos_fft(cand_f, cand_db, scale, unit_label)
if isempty(cand_f)
    return;
end
plot(cand_f, cand_db, 'wo', 'MarkerSize', 5, 'LineWidth', 1.0);
for i = 1:numel(cand_f)
    text(cand_f(i), cand_db(i), sprintf(' %.1f %s', cand_f(i) * scale, unit_label), ...
        'Color', [0.95 0.95 0.95], 'FontSize', 8, ...
        'VerticalAlignment', 'bottom');
end
end

function [cand_f, cand_db] = top_spectral_candidates(f, P_db, mask, n_cands, min_sep_hz)
idx_band = find(mask(:).' & isfinite(P_db(:).'));
idx_peaks = [];
for ii = 1:numel(idx_band)
    k = idx_band(ii);
    if k > 1 && k < numel(P_db) && P_db(k) >= P_db(k-1) && P_db(k) >= P_db(k+1)
        idx_peaks(end+1) = k; %#ok<AGROW>
    end
end
if isempty(idx_peaks)
    idx_peaks = idx_band;
end

[~, order] = sort(P_db(idx_peaks), 'descend');
idx_sorted = idx_peaks(order);
chosen = [];
for i = 1:numel(idx_sorted)
    k = idx_sorted(i);
    if isempty(chosen) || all(abs(f(k) - f(chosen)) >= min_sep_hz)
        chosen(end+1) = k; %#ok<AGROW>
    end
    if numel(chosen) >= n_cands
        break;
    end
end

cand_f = f(chosen);
cand_db = P_db(chosen);
end

function plot_fft_complex(x, Fs, nfft_max, xlim_values, unit_scale)
x = x(:);
N = min(numel(x), nfft_max);
Nfft = 2^nextpow2(N);
w = hamming(N);
X = fftshift(fft(x(1:N) .* w, Nfft));
f = (-Nfft/2:Nfft/2-1) * Fs / Nfft;
plot(f / unit_scale, 20*log10(abs(X) + 1e-12), 'LineWidth', 1.0);
xlim(xlim_values);
end

function [f, P_db] = single_sided_fft_db(x, Fs, pad_mult)
x = x(:);
N = numel(x);
Nfft = 2^nextpow2(max(1024, pad_mult * N));
X = fft((x - mean(x)) .* hamming(N), Nfft);
K = floor(Nfft/2) + 1;
f = (0:K-1) * Fs / Nfft;
P_db = 20*log10(abs(X(1:K)) + 1e-12);
end
