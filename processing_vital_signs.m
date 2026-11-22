function [rpm_est, bpm_est, info] = processing_vital_signs(buffer_IQ, Fs, dHP, dLP_rr)
% =========================================================================
% ECHO - Experimental Continuous-Wave Human Observation
% Copyright (C) 2026 Jose A. Hancco Mamani
% SPDX-License-Identifier: GPL-3.0-only
%
% This program is free software: you can redistribute it
% and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation,
% version 3 of the License.
% =========================================================================

% PROCESSING_VITAL_SIGNS Estimate respiratory and heart rates
% from complex baseband I/Q signals.
%
% Syntax:
%   [rpm_est, bpm_est, info] = processing_vital_signs( ...
%       buffer_IQ, Fs, dHP, dLP_rr)
%
% Description:
%   This function processes complex baseband I/Q samples
%   acquired by the ECHO SDR platform to estimate respiratory
%   rate (RR) and heart rate (HR).
%
%   The processing pipeline includes residual Tx-Rx carrier
%   offset compensation, dynamic I/Q circle centering,
%   phase extraction, cubic detrending, digital filtering,
%   spectral analysis, and confidence-based estimation.
%
% Inputs:
%   buffer_IQ - Complex baseband I/Q samples. The DC/clutter
%               component is retained for I/Q centering.
%   Fs        - Sampling frequency in Hz.
%   dHP       - High-pass digital filter object applied to
%               the extracted phase signal.
%   dLP_rr    - Low-pass digital filter object used to obtain
%               the respiratory signal.
%               Both filters must match Fs.
%
% Outputs:
%   rpm_est - Estimated respiratory rate in breaths/min.
%             NaN indicates that no reliable estimate
%             passed the confidence criteria.
%   bpm_est - Estimated heart rate in beats/min.
%             NaN indicates that no reliable estimate
%             passed the confidence criteria.
%   info    - Structure containing intermediate signals,
%             spectral estimates, and quality metrics used
%             for processing and diagnostic evaluation.
%
% Processing stages:
%   1. Input and filter validation.
%   2. Residual Tx-Rx carrier-offset compensation.
%   3. I/Q circle centering and phase extraction.
%   4. Cubic detrending, filtering, and edge trimming.
%   5. Respiratory-rate estimation.
%   6. Heart-rate estimation and confidence assessment.
%
% Requirements:
%   MATLAB and Signal Processing Toolbox.
%
% Note:
%   The function returns NaN when an estimate does not
%   satisfy the implemented confidence criteria.

% -------------------------------------------------------------------------
% 1. Analysis Parameters
% -------------------------------------------------------------------------
pad_mult    = 8;    % FFT zero-padding factor.
trim_s      = 5;    % Initial edge-trimming duration (s).
circ_win_s  = 15;   % Sliding-window duration for dynamic I/Q circle fitting (s).
hr_win_s    = 20;   % Window duration for HR consistency analysis (s).
welch_win_s = 15;   % Window duration for median-aggregated HR spectrum (s).

% Physiological search ranges (Hz).
rr_min_hz = 6/60;       % Minimum respiratory rate: 6 breaths/min.
rr_max_hz = 45/60;      % Maximum respiratory rate: 45 breaths/min.

hr_min_hz = 55/60;      % Minimum heart rate: 55 beats/min.
hr_max_hz = 160/60;     % Maximum heart rate: 160 beats/min.

% -------------------------------------------------------------------------
% 2. Input and Filter Validation
% -------------------------------------------------------------------------
if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('Fs must be a positive scalar in Hz.');
end

buffer_IQ = buffer_IQ(:);
if isempty(buffer_IQ)
    error('buffer_IQ is empty.');
end

if any(~isfinite(buffer_IQ))
    idx_good = isfinite(buffer_IQ);
    if nnz(idx_good) < 2
        error('buffer_IQ does not have enough finite samples.');
    end
    idx_all = (1:numel(buffer_IQ)).';
    buffer_IQ = interp1(idx_all(idx_good), buffer_IQ(idx_good), ...
        idx_all, 'linear', 'extrap');
end

check_filter_fs(dHP, 'dHP', Fs);
check_filter_fs(dLP_rr, 'dLP_rr', Fs);

dur_cap_s = numel(buffer_IQ)/Fs;
if dur_cap_s < 15
    warning('Captured in %.1f s: RR/HR may be somewhat unstable.', dur_cap_s);
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
    'rr_candidate_rpm', NaN, ...
    'rr_window_rpm', [], ...
    'rr_window_prominence', [], ...
    'rr_window_agreement', NaN, ...
    'rr_amplitude_ratio', NaN, ...
    'rr_consensus_prominence', NaN, ...
    'rr_signal_rms_rad', NaN, ...
    'hr_fallback', false, ...
    'hr_low_confidence', false, ...
    'hr_prominence', NaN, ...
    'hr_score_ratio', NaN, ...
    'hr_candidate_bpm', NaN, ...
    'hr_spectral_bpm_raw', NaN, ...
    'hr_spectral_bpm_fundamental', NaN, ...
    'hr_subharmonic_ratio', NaN, ...
    'hr_autocorr_bpm', NaN, ...
    'hr_autocorr_peak', NaN, ...
    'hr_folded_bpm', NaN, ...
    'hr_folded_agreement', NaN, ...
    'hr_folded_spread', NaN, ...
    'hr_signal_rms_rad', NaN, ...
    'hr_fusion_count', 0, ...
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
% 3. Residual Tx–Rx Carrier-Frequency Offset Compensation
% ========================================================================
S_iq = buffer_IQ;
if dur_cap_s >= 30
    [S_iq, f_off, vueltas, deriva_pp] = quitar_offset_portadora(S_iq, Fs);
else
    [S_iq, f_off, vueltas, deriva_pp] = quitar_offset_portadora(S_iq, Fs, 0);
end
info.carrier_offset_hz = f_off;
info.lo_drift_pp_rad = deriva_pp;

fprintf('\n--- RESIDUAL CARRIER OFFSET ---\n');
fprintf('Estimated TX-RX Offset       : %+.4f Hz\n', f_off);
fprintf('Total phase rotations          : %.2f over %.1f s\n', vueltas, dur_cap_s);
fprintf('Removed slow drift       : %.2f rad peak-to-peak\n', deriva_pp);

% ========================================================================
% 4. I/Q Circle Centering and Phase Extraction
% ========================================================================
[c_t, c0, R0, n_fits] = estimar_centro_circulo(S_iq, Fs, circ_win_s);
info.circle_radius = R0;
info.center_offset_ratio = abs(c0)/(R0 + eps);

Z = S_iq - c_t;
phi_raw = unwrap(angle(Z));
info.arc_span_rad = max(phi_raw) - min(phi_raw);

fprintf('\n--- I/Q GEOMETRY ---\n');
fprintf('Radius of the arc R            : %.5f\n', R0);
fprintf('|center| / R                : %.1f\n', info.center_offset_ratio);
fprintf('Angular opening of the arch   : %.2f rad\n', info.arc_span_rad);
fprintf('Valid circle windows : %d\n', n_fits);

if info.center_offset_ratio > 100
    warning(['|c|/R = %.0f: Clutter/leak dominates the utility signal. Check ' ...
        'antenna isolation, geometry, and orientation.'], info.center_offset_ratio);
end
if info.arc_span_rad < 0.25
    warning('IQ %.2f rad: very slight thoracic modulation.', info.arc_span_rad);
end
if info.arc_span_rad > 6
    warning('IQ opening of %.2f rad: likely to result in significant body movement.', info.arc_span_rad);
end

n = (0:numel(phi_raw)-1).';
radio = abs(Z);
radio_min = 0.20 * median(radio);
idx_ok = isfinite(phi_raw) & radio > radio_min;
if nnz(idx_ok) >= max(16, round(0.5*numel(phi_raw))) && ~all(idx_ok)
    phi_raw = interp1(n(idx_ok), phi_raw(idx_ok), n, 'linear', 'extrap');
end

% ========================================================================
% 5. Cubic Detrending, Filtering, and Edge Trimming
% ========================================================================
% Normalizing the sample-index axis improves numerical
% conditioning in the cubic polynomial fit.
u = linspace(-1, 1, numel(phi_raw)).';
p = polyfit(u, phi_raw, 3);
phi_det = phi_raw - polyval(p, u);


x_hp_full = filtfilt(dHP, phi_det);
x_rr_full = filtfilt(dLP_rr, x_hp_full);

trim_s = max(1.5, min(trim_s, 0.12 * dur_cap_s));
n_trim = round(trim_s * Fs);
if numel(x_hp_full) > 2*n_trim + round(10*Fs)
    keep = (n_trim+1):(numel(x_hp_full)-n_trim);
else
    warning('The capture is too short to crop %.1f s from each edge.', trim_s);
    keep = 1:numel(x_hp_full);
end

x = x_hp_full(keep);
x_rr = x_rr_full(keep);
phi_use = phi_raw(keep);

info.phi = phi_use(:);
info.x = x(:);
info.x_rr = x_rr(:);
info.rr_signal_rms_rad = sqrt(mean((x_rr-median(x_rr)).^2));

fprintf('Signal: %.1f s (%d samples at %.1f Hz)\n', numel(x)/Fs, numel(x), Fs);

% ========================================================================
% 6. Respiratory-Rate Estimation
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
    f_rr_global = f_rr_est;
    info.f_rr_est = f_rr_est;
    info.rr_candidate_rpm = double(60*f_rr_est);

    [info.rr_window_rpm, info.rr_window_prominence, ...
        info.rr_window_agreement, info.rr_amplitude_ratio, ...
        f_rr_consensus, info.rr_consensus_prominence] = ...
        consistencia_rr(x_rr, Fs, f_rr_est, rr_min_hz, rr_max_hz, pad_mult);

    n_rr_win = numel(info.rr_window_rpm);
    if isfinite(f_rr_consensus)
        info.f_rr_est = f_rr_consensus;
        info.rr_candidate_rpm = double(60*f_rr_consensus);
    end

    evidencia_global = rr_quality.prominence >= 2.0 && ...
                       rr_quality.score_ratio >= 1.10 && ...
                       abs(60*f_rr_global-info.rr_candidate_rpm) <= ...
                           max(3.5, 0.18*info.rr_candidate_rpm);
    evidencia_ventanas = info.rr_consensus_prominence >= 2.0;
    acuerdo_rr_requerido = 0.75;
    if n_rr_win == 3
        acuerdo_rr_requerido = 2/3;
    end
    acuerdo_rr_ok = info.rr_window_agreement+eps >= acuerdo_rr_requerido;
   
    evidencia_global_fuerte = rr_quality.prominence >= 3.0 && ...
        rr_quality.score_ratio >= 1.50;
    evidencia_rr_lenta = info.rr_candidate_rpm <= 8 && ...
        rr_quality.prominence >= 3.0 && rr_quality.score_ratio >= 1.15 && ...
        info.rr_amplitude_ratio >= 0.30;
    decision_espectral_rr = (acuerdo_rr_ok && ...
        (evidencia_global || evidencia_ventanas)) || ...
        evidencia_global_fuerte || evidencia_rr_lenta;

    info.rr_low_confidence = n_rr_win < 3 || ...
                             info.rr_amplitude_ratio < 0.15 || ...
                             info.rr_signal_rms_rad < 0.15 || ...
                             ~decision_espectral_rr;

    fprintf('\n--- RR CONSISTENCY ---\n');
    fprintf('RPM per window (rpm): %s\n', formatear_lista(info.rr_window_rpm));
    fprintf('RR Prominence per Window: %s\n', ...
        formatear_lista(info.rr_window_prominence));
    fprintf('Agreement with a Candidate (+-tolerance): %.0f%%\n', ...
        100*info.rr_window_agreement);
    fprintf('Min/Max Amplitude Persistence: %.2f\n', ...
        info.rr_amplitude_ratio);
    fprintf('Median prominence (consensus): %.2f\n', ...
        info.rr_consensus_prominence);
    fprintf('Respiratory Phase RMS       : %.4f rad\n', ...
        info.rr_signal_rms_rad);
    if info.rr_low_confidence
        rpm_est = NaN;
        warning(['RR candidate %.1f rpm REJECTED due to low confidence ' ...
            '(prom %.2f, ratio %.2f, agreement %.0f%%, persistence %.2f).'], ...
            info.rr_candidate_rpm, rr_quality.prominence, ...
            rr_quality.score_ratio, 100*info.rr_window_agreement, ...
            info.rr_amplitude_ratio);
    else
        rpm_est = info.rr_candidate_rpm;
    end
end

% ========================================================================
% 7. Heart-Rate Estimation
% ========================================================================
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

rr_for_hr = NaN;

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

hr_bpm_raw = double(60*f_hr_est);
info.hr_spectral_bpm_raw = hr_bpm_raw;
info.hr_subharmonic_ratio = soporte_subarmonico( ...
    f_hr_ax, P_hr, f_hr_est/2, f_hr_est, 0.035);

hr_bpm_espectral = hr_bpm_raw;
if f_hr_est/2 >= hr_min_hz && info.hr_subharmonic_ratio >= 0.80
    hr_bpm_espectral = hr_bpm_raw/2;
end
info.hr_spectral_bpm_fundamental = hr_bpm_espectral;

[info.hr_window_bpm, ~] = consistencia_hr( ...
    x_hr, Fs, hr_win_s, hr_min_hz, hr_max_hz, hr_bpm_raw, rr_for_hr);
n_hr_win = numel(info.hr_window_bpm);

[info.hr_folded_bpm, info.hr_folded_agreement, ...
    info.hr_folded_spread] = consenso_hr_con_armonicos( ...
        info.hr_window_bpm, 60*hr_min_hz, 60*hr_max_hz);


[info.hr_autocorr_bpm, info.hr_autocorr_peak] = ...
    candidato_hr_autocorr(x_hr, Fs, 60*hr_min_hz, 60*hr_max_hz);

tol_fusion_bpm = 7;
est = [hr_bpm_espectral, info.hr_autocorr_bpm, info.hr_folded_bpm];
par_espectral_ac = cercanos(est(1), est(2), tol_fusion_bpm);
par_espectral_ventanas = cercanos(est(1), est(3), tol_fusion_bpm) && ...
    info.hr_folded_agreement >= 0.75 && info.hr_folded_spread <= 5;
par_ac_ventanas = cercanos(est(2), est(3), tol_fusion_bpm) && ...
    info.hr_folded_agreement >= 0.50 && info.hr_folded_spread <= 3 && ...
    info.hr_autocorr_peak >= 0.22;

usar = false(1,3);
if par_espectral_ac, usar([1 2]) = true; end
if par_espectral_ventanas, usar([1 3]) = true; end
if par_ac_ventanas, usar([2 3]) = true; end
usar = usar & isfinite(est);
info.hr_fusion_count = nnz(usar);
if any(usar)
    hr_bpm_fusion = median(est(usar));
else
    hr_bpm_fusion = NaN;
end

x_hr_c = x_hr - median(x_hr);
info.hr_signal_rms_rad = sqrt(mean(x_hr_c.^2));

presencia_hr = info.hr_signal_rms_rad >= 0.010;
periodicidad_hr = info.hr_autocorr_peak >= 0.20;
evidencia_espectral_hr = hr_quality.prominence >= 1.40;
fusion_valida = isfinite(hr_bpm_fusion) && info.hr_fusion_count >= 2;

info.hr_low_confidence = n_hr_win < 3 || ~presencia_hr || ...
    ~periodicidad_hr || ~evidencia_espectral_hr || ~fusion_valida;

if fusion_valida
    info.hr_candidate_bpm = hr_bpm_fusion;
    info.f_hr_est = hr_bpm_fusion/60;
    info.hr_agreement = acuerdo_hr_con_armonicos( ...
        info.hr_window_bpm, hr_bpm_fusion, 5);
else
    info.hr_candidate_bpm = hr_bpm_raw;
    info.f_hr_est = f_hr_est;
    info.hr_agreement = info.hr_folded_agreement;
end

fprintf('\n--- HR CONSISTENCY ---\n');
fprintf('HR per window (bpm): %s\n', formatear_lista(info.hr_window_bpm));
fprintf('Raw/Fundamental Spectral HR: %.1f / %.1f bpm\n', ...
    info.hr_spectral_bpm_raw, info.hr_spectral_bpm_fundamental);
fprintf('f/2 Subharmonic Support      : %.2f\n', info.hr_subharmonic_ratio);
fprintf('HR by autocorrelation       : %.1f bpm (peak %.2f)\n', ...
    info.hr_autocorr_bpm, info.hr_autocorr_peak);
fprintf('Temporary consensus reached    : %.1f bpm (agreement %.0f%%, dispersion %.1f)\n', ...
    info.hr_folded_bpm, 100*info.hr_folded_agreement, info.hr_folded_spread);
fprintf('Cardiac Phase RMS         : %.4f rad\n', info.hr_signal_rms_rad);
if fusion_valida
    fprintf('Fusion HR (%d estimators)    : %.1f bpm\n', ...
        info.hr_fusion_count, hr_bpm_fusion);
else
    fprintf('Fusion HR                    : no consensus\n');
end
if isfinite(info.hr_agreement)
    fprintf('Harmonic-aware agreement with the Candidate (+-5 bpm): %.0f%%\n', ...
        100*info.hr_agreement);
else
    fprintf('Agreement with candidate: cannot be evaluated (very few windows)\n');
end

if info.hr_low_confidence
    bpm_est = NaN;
    warning(['HR candidate %.1f bpm REJECTED due to low confidence ' ...
        '(prom %.2f, RMS %.4f rad, AC %.2f, fusion %d).'], ...
        info.hr_candidate_bpm, hr_quality.prominence, ...
        info.hr_signal_rms_rad, info.hr_autocorr_peak, ...
        info.hr_fusion_count);
else
    bpm_est = hr_bpm_fusion;
end
end

% =========================================================================
% Auxiliary Functions
% =========================================================================
% Window consistency algorithms, Kasa circle fitting, and robust peak selection ---
% ------------------------------------------------------------------------
% Respiratory-Rate Consistency Assessment
% Evaluates RR estimates across overlapping signal windows.
% Computes spectral prominence, temporal agreement, and
% amplitude persistence to assess estimation reliability.
% ------------------------------------------------------------------------
function [rpm_win, prom_win, acuerdo, amp_ratio, f_consenso, prom_consenso] = consistencia_rr( ...
    x_rr, Fs, f_ref, fmin, fmax, pad_mult)
x_rr = x_rr(:);
N = numel(x_rr);
dur_s = N/Fs;
win_s = min(20, max(12, 0.40*dur_s));
Nw = round(win_s*Fs);
hop = max(1, round(Nw/2));
rpm_win = [];
prom_win = [];
amp_ref = [];

if Nw >= N
    acuerdo = 0;
    amp_ratio = 0;
    f_consenso = NaN;
    prom_consenso = 0;
    return;
end

for i0 = 1:hop:(N-Nw+1)
    seg = x_rr(i0:i0+Nw-1);
    [f, P, bin_hz] = espectro(seg, Fs, pad_mult);
    mask = (f >= fmin) & (f <= fmax);
    [fw, qw] = select_rr_peak_robust(f, P, mask, bin_hz);
    rpm_win(end+1,1) = 60*fw;
    prom_win(end+1,1) = qw.prominence;

    t = (0:Nw-1).'/Fs;
    A = [ones(Nw,1), t, sin(2*pi*f_ref*t), cos(2*pi*f_ref*t)];
    q = A\seg;
    amp_ref(end+1,1) = hypot(q(3),q(4));
end

validas = isfinite(rpm_win) & prom_win >= 1.2;

if any(validas)
    rpm_mediana = median(rpm_win(validas));
    tol_rpm = max(3.5, 0.18*rpm_mediana);
    mejor = validas & abs(rpm_win-rpm_mediana) <= tol_rpm;
else
    mejor = false(size(rpm_win));
end

acuerdo = nnz(mejor)/max(numel(rpm_win),1);
if any(mejor)
    rpm_consenso = median(rpm_win(mejor));
    f_consenso = rpm_consenso/60;
    prom_consenso = median(prom_win(mejor));
else
    f_consenso = NaN;
    prom_consenso = 0;
end

amp_ref = amp_ref(isfinite(amp_ref) & amp_ref > 0);
if isempty(amp_ref)
    amp_ratio = 0;
else
    amp_ratio = min(amp_ref)/max(amp_ref);
end
end

% =========================================================================
% I/Q Circle Fitting and Dynamic Center Estimation
% =========================================================================
% ------------------------------------------------------------------------
% Dynamic I/Q Circle Center Estimation
% Estimates the global I/Q circle and tracks its center
% across overlapping windows using Kasa circle fitting.
% Valid window estimates are smoothed and interpolated
% to obtain a time-varying center trajectory.
% ------------------------------------------------------------------------
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

    centros(end+1,1) = cw;
    tiempos(end+1,1) = i0 + (Nw-1)/2;
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

% ------------------------------------------------------------------------
% Kasa Algebraic Circle Fitting
% Estimates the center and radius of a circle from
% complex I/Q samples using algebraic least squares.
% ------------------------------------------------------------------------
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
% Window-Based Heart-Rate Consistency
% =========================================================================
% ------------------------------------------------------------------------
% Heart-Rate Consistency Assessment
% Evaluates HR estimates across overlapping signal windows
% and measures their agreement with the reference candidate.
% ------------------------------------------------------------------------
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


    [f0, ~] = select_hr_peak_robust(f, P, m, f_rr);
    if isfinite(f0)
        bpm_win(end+1,1) = 60*f0;
    end
end

if numel(bpm_win) < 2
    acuerdo = NaN;
else
    acuerdo = mean(abs(bpm_win-bpm_ref) <= 5);
end
end

function ratio = soporte_subarmonico(f, P, f_sub, f_full, half_width)
if ~isfinite(f_sub) || ~isfinite(f_full) || f_sub <= 0 || f_full <= 0
    ratio = 0;
    return;
end
m_sub = abs(f-f_sub) <= half_width;
m_full = abs(f-f_full) <= half_width;
if ~any(m_sub) || ~any(m_full)
    ratio = 0;
    return;
end
ratio = max(P(m_sub))/(max(P(m_full))+eps);
end

function [bpm, pico] = candidato_hr_autocorr(x, Fs, bpm_min, bpm_max)
x = x(:) - mean(x);
if numel(x) < round(3*Fs) || sqrt(mean(x.^2)) <= eps
    bpm = NaN;
    pico = 0;
    return;
end
r = xcorr(x, 'coeff');
r = r(numel(x):end);
i_lo = max(2, round(60/bpm_max*Fs)+1);
i_hi = min(numel(r)-1, round(60/bpm_min*Fs)+1);
if i_hi <= i_lo
    bpm = NaN;
    pico = 0;
    return;
end
idx = i_lo:i_hi;
loc = idx(r(idx) >= r(idx-1) & r(idx) >= r(idx+1));
if isempty(loc)
    [pico, k] = max(r(idx));
    kk = idx(k);
else
    [pico, k] = max(r(loc));
    kk = loc(k);
end

d = 0;
if kk > 1 && kk < numel(r)
    den = r(kk-1) - 2*r(kk) + r(kk+1);
    if abs(den) > 1e-12
        d = 0.5*(r(kk-1)-r(kk+1))/den;
        d = max(-0.5,min(0.5,d));
    end
end
lag = (kk-1)+d;
bpm = 60*Fs/lag;
end

% ------------------------------------------------------------------------
% Harmonic-Aware Temporal Consensus
% Evaluates agreement among window-based HR estimates,
% considering both the detected frequency and its
% half-frequency alternative to account for harmonic ambiguity.
% ------------------------------------------------------------------------
function [centro, acuerdo, dispersion] = consenso_hr_con_armonicos( ...
    bpm_win, bpm_min, bpm_max)
bpm_win = bpm_win(isfinite(bpm_win));
if isempty(bpm_win)
    centro = NaN;
    acuerdo = NaN;
    dispersion = Inf;
    return;
end

opciones = [];
for k = 1:numel(bpm_win)
    if bpm_win(k) >= bpm_min && bpm_win(k) <= bpm_max
        opciones(end+1,:) = [bpm_win(k), k];
    end
    mitad = bpm_win(k)/2;
    if mitad >= bpm_min && mitad <= bpm_max
        opciones(end+1,:) = [mitad, k];
    end
end
if isempty(opciones)
    centro = NaN;
    acuerdo = 0;
    dispersion = Inf;
    return;
end

rejilla = bpm_min:0.25:bpm_max;
cuentas = zeros(size(rejilla));
for j = 1:numel(rejilla)
    votos = false(1,numel(bpm_win));
    for q = 1:size(opciones,1)
        if abs(opciones(q,1)-rejilla(j)) <= 5
            votos(opciones(q,2)) = true;
        end
    end
    cuentas(j) = nnz(votos);
end
[n_votos, j] = max(cuentas);
centro_previo = rejilla(j);

elegidos = NaN(size(bpm_win));
for k = 1:numel(bpm_win)
    vals = bpm_win(k);
    mitad = bpm_win(k)/2;
    if mitad >= bpm_min && mitad <= bpm_max
        vals = [vals, mitad];
    end
    [dist, q] = min(abs(vals-centro_previo));
    if dist <= 5
        elegidos(k) = vals(q);
    end
end
elegidos = elegidos(isfinite(elegidos));
acuerdo = n_votos/max(numel(bpm_win),1);
if isempty(elegidos)
    centro = NaN;
    dispersion = Inf;
else
    centro = median(elegidos);
    dispersion = 1.4826*median(abs(elegidos-centro));
end
end

function acuerdo = acuerdo_hr_con_armonicos(bpm_win, bpm_ref, tolerancia)
bpm_win = bpm_win(isfinite(bpm_win));
if isempty(bpm_win) || ~isfinite(bpm_ref)
    acuerdo = NaN;
    return;
end
ok = false(size(bpm_win));
for k = 1:numel(bpm_win)
    ok(k) = abs(bpm_win(k)-bpm_ref) <= tolerancia || ...
        abs(bpm_win(k)/2-bpm_ref) <= tolerancia;
end
acuerdo = mean(ok);
end

function tf = cercanos(a, b, tolerancia)
tf = isfinite(a) && isfinite(b) && abs(a-b) <= tolerancia;
end

% =========================================================================
% Spectrum Estimation
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
    Pmat(:,end+1) = Pi(:);
end
P = median(Pmat, 2);
end

% =========================================================================
% Robust Spectral Peak Selection
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

    scores(i) = prom;
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
if ~isfinite(f_est) || f_est < 10/60 || ...
        (f_est >= 12/60 && f_est < 24/60)
    return;
end

k_est = nearest_masked_bin(f_ax, mask, f_est);
if isempty(k_est), return; end
p_est = P(k_est);
best_f = f_est;
best_ratio = 0;

for divisor = [2 3]
    f_candidate = f_est/divisor;
    if f_candidate < 5.5/60 || f_candidate > 24/60
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

    scores(i) = prom/(1+harm_pen);
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
        idx_peaks(end+1) = k;
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
    txt = 'no data';
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
    error('%s was designed for %.3f Hz, but the current Fs is %.3f Hz.', ...
        name, fs_design, fs_runtime);
end
end

% =========================================================================
% Residual Carrier-Offset Estimation and Compensation
% =========================================================================
% ------------------------------------------------------------------------
% Residual Carrier-Frequency Offset Compensation
% Estimates the residual frequency offset using FFT analysis
% and phase-slope refinement, then compensates the resulting
% phase rotation. Optional low-frequency phase-drift removal
% is controlled by the fc_track parameter.
% ------------------------------------------------------------------------
function [z_out, f_off, vueltas, deriva_pp] = quitar_offset_portadora(z, Fs, fc_track)
if nargin < 3 || isempty(fc_track), fc_track = 0; end

N = numel(z);
n = (0:N-1).';

Nf = 2^nextpow2(8*N);
Zf = fft(z, Nf);
f = (0:Nf-1).'/Nf*Fs;
f(f > Fs/2) = f(f > Fs/2)-Fs;
mag = abs(Zf);
mag(abs(f) > 2) = 0;
[~, k] = max(mag);
f0 = f(k);

z1 = z .* exp(-1j*2*pi*f0*n/Fs);

psi = unwrap(angle(z1));
u = linspace(-1,1,N).';
p = polyfit(u, psi, 1);
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
