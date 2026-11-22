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

% =========================================================================
% ECHO - Main Processing and Visualization Script
%
% Description:
%   Executes the main ECHO processing pipeline, from raw I/Q data
%   acquisition to physiological parameter estimation and visualization.
%
%   The script coordinates the following modules:
%       1. cw.m - Continuous-wave radar data acquisition.
%       2. frontend_dsp.m - Low-IF downconversion and decimation.
%       3. processing_vital_signs.m - Respiratory and heart rate estimation.
%
% Requirements:
%   - MATLAB
%   - Signal Processing Toolbox
%   - bladeRF MATLAB interface
%   - Nuand bladeRF 2.0 micro xA4 SDR
%
% Usage:
%   Run this script after placing the required ECHO modules
%   in the MATLAB working directory.
%
% =========================================================================
% Initialize the workspace by clearing variables, the command window, and closing any open figures
clear; clc; close all;

% =========================================================================
% 1. I/Q DATA ACQUISITION
% Execute the CW acquisition module to obtain raw complex I/Q
% samples and the associated hardware configuration metadata.
% =========================================================================
try
    [rx_continuous, cfg] = cw();
catch ME
    error('Error in cw.m: %s', ME.message);
end

if isempty(rx_continuous)
    error('cw.m does not generate samples.');
end
% Abort processing if streaming overruns are detected, as dropped
% samples can compromise phase continuity.
if isfield(cfg, 'overruns') && cfg.overruns > 0
    error(['The capture had %d overruns. It is not being processed because samples are missing. ' ...
        'and the phase ceases to be continuous.'], cfg.overruns);
end

% =========================================================================
% 2. LOW-IF TO BASEBAND CONVERSION
% Apply digital downconversion and multi-stage decimation
% to obtain the complex baseband signal at the target
% sampling frequency.
% =========================================================================
% Downconvert the raw Low-IF signal to baseband and decimate it to a lower sampling rate (e.g., 100 Hz)
fprintf('\nRunning a front-end DSP...\n');
[iq_frontend, Fs_frontend] = frontend_dsp(rx_continuous, cfg);
if isempty(iq_frontend) || ~isfinite(Fs_frontend) || Fs_frontend <= 0
    error('frontend_dsp returned invalid output.');
end

% =========================================================================
% 3. PHYSIOLOGICAL FILTER DESIGN
% Check the available signal duration and design the high-pass
% and low-pass filters required for physiological signal analysis.
% =========================================================================
fc_hp = 0.04; % High-pass filter cutoff frequency (Hz).

dur_s = numel(iq_frontend)/Fs_frontend;
if dur_s < 25
    warning(['Capture at only %.1f s with HP %.2f Hz. For more stable estimates of ' ...
        'RR/HR, approximately 30 seconds is recommended.'], dur_s, fc_hp);
end

fprintf('Physiological filter: HP %.2f Hz, duration %.1f s.\n', fc_hp, dur_s);

% Design a high-pass IIR filter to suppress slow baseline drift
% and low-frequency stationary clutter.
dHP = designfilt('highpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', fc_hp, 'SampleRate', Fs_frontend);

% Design a low-pass IIR filter to limit the respiratory
% signal bandwidth.
dLP_rr = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 0.8, 'SampleRate', Fs_frontend);

% =========================================================================
% 4. PHYSIOLOGICAL PARAMETER ESTIMATION
% Process the baseband signal using the designed filters
% and estimate the respiratory rate (RR) and heart rate (HR).
% =========================================================================
[rpm, bpm, vital_info] = processing_vital_signs( ...
    iq_frontend, Fs_frontend, dHP, dLP_rr);


% Display the estimated physiological parameters and their
% corresponding confidence or rejection status.
fprintf('\n=============================\n');
fprintf('RESULTS\n');
if vital_info.rr_fallback || ~isfinite(rpm)
    if isfield(vital_info,'rr_candidate_rpm') && isfinite(vital_info.rr_candidate_rpm)
        fprintf('Respiration: NOT DETECTED (candidate %.1f rpm rejected)\n', vital_info.rr_candidate_rpm);
    else
        fprintf('Respiration: NOT DETECTED\n');
    end
elseif vital_info.rr_low_confidence
    fprintf('Respiration: %.1f rpm  [LOW CONFIDENCE]\n', rpm);
else
    fprintf('Respiration: %.1f rpm\n', rpm);
end

if vital_info.hr_fallback || ~isfinite(bpm)
    if isfield(vital_info,'hr_candidate_bpm') && isfinite(vital_info.hr_candidate_bpm)
        fprintf('Heart rate: NOT DETECTED (candidate %.1f bpm rejected)\n', vital_info.hr_candidate_bpm);
    else
        fprintf('Heart rate: NOT DETECTED\n');
    end
elseif vital_info.hr_low_confidence
    fprintf('Heart Rate: %.1f bpm  [LOW CONFIDENCE]\n', bpm);
else
    fprintf('Heart Rate: %.1f bpm\n', bpm);
end

if isfinite(vital_info.hr_agreement)
    fprintf('HR Agreement Between Windows: %.0f%%\n', 100*vital_info.hr_agreement);
end
fprintf('=============================\n');
fprintf('Current Final HR: %.4f bpm\n', bpm);
fprintf('Current Final RR: %.4f rpm\n', rpm);

% =========================================================================
% 5. DIAGNOSTIC VISUALIZATION
% Generate a diagnostic dashboard showing acquisition quality,
% baseband I/Q behavior, phase signals, and physiological estimates.
% =========================================================================
% Launch the diagnostic dashboard to visualize hardware metrics, phase extraction, and spectra
plot_diagnostico_radar(rx_continuous, cfg, iq_frontend, Fs_frontend, ...
    rpm, bpm, vital_info);

% -------------------------------------------------------------------------
% PLOT_DIAGNOSTICO_RADAR
% Generates a 3-by-3 diagnostic dashboard displaying raw signal
% quality, baseband I/Q behavior, phase signals, respiratory
% and cardiac spectra, and window-level HR estimates.
% -------------------------------------------------------------------------
function plot_diagnostico_radar(rx_raw, cfg, iq_bb, Fs_bb, rpm, bpm, info)
Fs_raw = cfg.fs;
rx_raw = rx_raw(:);
iq_bb = iq_bb(:);

rr_txt = 'NO DET.';
if isfinite(rpm), rr_txt = sprintf('%.1f rpm', rpm); end
hr_txt = 'NO DET.';
if isfinite(bpm), hr_txt = sprintf('%.1f bpm', bpm); end

figure('Name', 'Radar Vital Signs Diagnosis', ...
    'Units', 'normalized', 'Position', [0.04 0.06 0.92 0.86]);
sgtitle(sprintf('Radar | RR %s | HR %s | |c|/R %.2f | arc %.2f rad', ...
    rr_txt, hr_txt, info.center_offset_ratio, info.arc_span_rad), ...
    'FontWeight', 'bold', 'FontSize', 13);

% 1) RX level per block: Receiver amplitude and RMS levels per block to check for ADC saturation or missing data
subplot(3,3,1);
if isfield(cfg,'block_peak') && isfield(cfg,'block_rms') && ...
   isfield(cfg,'actual_counts') && isfield(cfg,'block_len')
    n_blk = numel(cfg.block_peak);
    t_blk = (0:n_blk-1)*(cfg.block_len/cfg.fs);
    plot(t_blk,cfg.block_peak,'b'); hold on;
    plot(t_blk,cfg.block_rms,'c');
    yline(0.99,':r');
    idx_inc = cfg.actual_counts ~= cfg.block_len;
    if any(idx_inc)
        plot(t_blk(idx_inc),cfg.block_peak(idx_inc),'rx');
    end
    grid on; xlabel('s'); ylabel('amplitude');
    title(sprintf('RX (%d incomplete blocks)', nnz(idx_inc)));
else
    n_raw = min(numel(rx_raw), round(0.01*Fs_raw));
    t_ms = (0:n_raw-1)/Fs_raw*1e3;
    plot(t_ms,real(rx_raw(1:n_raw))); hold on;
    plot(t_ms,imag(rx_raw(1:n_raw)));
    grid on; xlabel('ms'); ylabel('amplitude'); title('Raw RX');
end

% 2) IQ without centering: Scatter plot of the baseband I/Q constellation to visualize the cardiopulmonary phase arc
subplot(3,3,2);
idx_c = round(linspace(1,numel(iq_bb),min(8000,numel(iq_bb))));
plot(real(iq_bb(idx_c)),imag(iq_bb(idx_c)),'.','MarkerSize',3); hold on;
plot(0,0,'r+','MarkerSize',10,'LineWidth',1.5);
axis equal; grid on; xlabel('I'); ylabel('Q');
title(sprintf('Base-band IQ @ %.0f Hz', Fs_bb));

% 3) Raw low-IF FFT
subplot(3,3,3);
plot_fft_complex(rx_raw,Fs_raw,262144,[-150 150],1e3);
xline(cfg.tx_tone_hz/1e3,'--r',sprintf('%.0f kHz',cfg.tx_tone_hz/1e3));
grid on; xlabel('kHz'); ylabel('dB'); title('Raw FFT RX');

% 4) Phase
subplot(3,3,4);
t = (0:numel(info.phi)-1).'/info.Fs;
plot(t,info.phi,'LineWidth',0.8);
grid on; xlabel('s'); ylabel('rad'); title('Unwrapped phase');

% 5) Respiration
subplot(3,3,5);
t = (0:numel(info.x_rr)-1).'/info.Fs;
plot(t,info.x_rr,'LineWidth',1.0);
grid on; xlabel('s'); ylabel('rad'); title(['Respiration: ' rr_txt]);

% 6) Heart rate
subplot(3,3,6);
t = (0:numel(info.x_hr)-1).'/info.Fs;
plot(t,info.x_hr,'LineWidth',0.8);
grid on; xlabel('s'); ylabel('rad'); title(['Cardiac: ' hr_txt]);

% 7) Spectrum RR: Frequency spectra for respiration and heart rate, highlighting the selected peaks and harmonics
subplot(3,3,7);
plot(info.f_rr_axis,20*log10(info.P_rr+1e-12),'LineWidth',1.1); hold on;
if isfinite(info.f_rr_est)
    xline(info.f_rr_est,'--r',rr_txt);
end
grid on; xlim([0 1.2]); xlabel('Hz'); ylabel('dB'); title('Respiratory spectrum');

% 8) Spectrum HR: Frequency spectra for respiration and heart rate, highlighting the selected peaks and harmonics
subplot(3,3,8);
plot(info.f_hr_axis,20*log10(info.P_hr+1e-12),'LineWidth',1.1); hold on;
if isfinite(info.f_hr_est)
    xline(info.f_hr_est,'--r',hr_txt);
end
if isfinite(info.f_rr_est)
    for h = 2:6
        fh = h*info.f_rr_est;
        if fh >= 0.4 && fh <= 3.6
            xline(fh,':',sprintf('%dxRR',h),'Color',[0.6 0.6 0.6]);
        end
    end
end
grid on; xlim([0.4 3.6]); xlabel('Hz'); ylabel('dB'); title('Cardiac spectrum');

% 9) HR per windows
subplot(3,3,9);
if ~isempty(info.hr_window_bpm)
    stem(1:numel(info.hr_window_bpm),info.hr_window_bpm,'filled'); hold on;
    if isfinite(bpm), yline(bpm,'--r'); end
    grid on; xlabel('window'); ylabel('bpm');
    if isfinite(info.hr_agreement)
        title(sprintf('HR per window (agreement %.0f%%)',100*info.hr_agreement));
    else
        title('HR per window');
    end
else
    text(0.5,0.5,'No HR windows available','HorizontalAlignment','center'); axis off;
end
end

% -------------------------------------------------------------------------
% PLOT_FFT_COMPLEX
% Computes and displays the centered magnitude spectrum of
% a complex input signal using the FFT, with magnitude in dB.
% -------------------------------------------------------------------------
function plot_fft_complex(x, Fs, nfft_max, xlim_values, unit_scale)
x = x(:);
N = min(numel(x),nfft_max);
Nfft = 2^nextpow2(N);
w = hamming(N);
X = fftshift(fft(x(1:N).*w,Nfft));
f = (-Nfft/2:Nfft/2-1)*Fs/Nfft;
plot(f/unit_scale,20*log10(abs(X)+1e-12),'LineWidth',1.0);
xlim(xlim_values);
end
