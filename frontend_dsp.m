function [iq_out, Fs, front_info] = frontend_dsp(iq_in, cfg)
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

% FRONTEND_DSP Converts the received Low-IF I/Q signal to complex
% baseband and performs multi-stage decimation.
%
% Syntax:
%   [iq_out, Fs, front_info] = frontend_dsp(iq_in, cfg)
%
% Description:
%   This function processes complex I/Q samples acquired by the
%   ECHO SDR platform. It estimates the actual received Low-IF
%   tone, digitally downconverts the signal using the measured
%   frequency, and applies three FIR decimation stages.
%
% Inputs:
%   iq_in - Complex Low-IF I/Q input samples.
%   cfg   - Configuration structure containing:
%       fs         - Input sampling frequency (Hz).
%       tx_tone_hz - Nominal Low-IF tone frequency (Hz).
%       debug_plots - Optional flag to enable diagnostic plots.
%
% Outputs:
%   iq_out    - Complex baseband I/Q samples after filtering
%               and decimation.
%   Fs        - Output sampling frequency (Hz).
%   front_info - Structure containing Low-IF frequency
%                diagnostics:
%       nominal_hz       - Nominal tone frequency (Hz).
%       measured_hz      - Measured tone frequency (Hz).
%       offset_hz        - Frequency deviation (Hz).
%       contrast_db      - Peak-to-local-median spectral contrast (dB).
%       exceeded_2khz    - Flag indicating an offset strictly
%                          greater than 2 kHz.
%       mix_frequency_hz - Frequency used for digital mixing (Hz).
%
% Processing stages:
%   1. Input validation and preparation.
%   2. Low-IF tone estimation and spectral contrast assessment.
%   3. Digital downconversion using the measured tone frequency.
%   4. Multi-stage FIR filtering and decimation.
%
% Notes:
%   The measured tone frequency is used for digital mixing,
%   including when its deviation from nominal exceeds 2 kHz.
%   This condition generates a warning but does not reject
%   the acquisition.
%
%   The total decimation factor is 10,000. An input sampling
%   frequency of 1 MHz therefore produces a 100 Hz output.

% Validate the required configuration fields.
if nargin < 2 || ~isstruct(cfg) || ~isfield(cfg,'fs') || ~isfield(cfg,'tx_tone_hz')
    error('The cfg file must contain fs and tx_tone_hz.');
end

Fs = cfg.fs;

% Validate the input sampling frequency.
if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('cfg.fs must be positive.');
end
% Reshape the input I/Q samples as a column vector.
iq = iq_in(:);

% Initialize the Low-IF frequency diagnostic structure.
front_info = struct();
front_info.nominal_hz = cfg.tx_tone_hz;
front_info.measured_hz = NaN;
front_info.offset_hz = NaN;
front_info.contrast_db = NaN;
front_info.exceeded_2khz = false;
front_info.mix_frequency_hz = NaN;
% Check for empty input and handle non-finite samples.
if isempty(iq)
    error('iq_in is empty.');
end

if any(~isfinite(iq))
    good = isfinite(iq);
    if nnz(good) < 2
        error('iq_in must contain at least two finite samples.');
    end
    idx = (1:numel(iq)).';
    iq = interp1(idx(good),iq(good),idx,'linear','extrap');
end

debug_plots = isfield(cfg,'debug_plots') && logical(cfg.debug_plots);

% Remove the complex mean to suppress the DC component
% before estimating the Low-IF tone.
iq = iq - mean(iq);

if debug_plots
    plot_spectrum(iq,Fs,65536,'FFT before mixing','kHz',1e3);
end
% Estimate the received Low-IF tone and its spectral contrast.
[f_tone, tone_contrast_db] = estimar_tono(iq,Fs,cfg.tx_tone_hz);

% Store the measured frequency and associated diagnostics.
front_info.measured_hz = f_tone;
front_info.offset_hz = f_tone - cfg.tx_tone_hz;
front_info.contrast_db = tone_contrast_db;

% Flag frequency deviations strictly greater than 2 kHz.
% This is a diagnostic flag and does not reject the signal.
front_info.exceeded_2khz = ...
    abs(front_info.offset_hz) > 2000;

% Use the measured tone frequency for digital downconversion,
% even when the deviation exceeds 2 kHz.
front_info.mix_frequency_hz = f_tone;
fprintf('Frontend: nominal frequency %.1f Hz, measured %.2f Hz (offset %+.2f Hz)\n', ...
    cfg.tx_tone_hz,f_tone,f_tone-cfg.tx_tone_hz);
fprintf('Frontend: tone contrast of %.1f dB above the local median\n',tone_contrast_db);
if tone_contrast_db < 10
    warning(['The Low-IF tone has only %.1f dB of contrast. Check TX, ' ...
        'antennas, interference or saturation before relying on RR/HR.'], ...
        tone_contrast_db);
end
% Generate the sample-index vector for the digital mixing oscillator.
n = (0:numel(iq)-1).'; 

% Shift the measured Low-IF tone to complex baseband.
iq = iq .* exp(-1j*2*pi*f_tone*n/Fs);
% -------------------------------------------------------------------------
% Multi-stage FIR decimation.
% For an input sampling rate of 1 MHz:
% 1 MHz -> 20 kHz -> 1 kHz -> 100 Hz.
% -------------------------------------------------------------------------

% Stage 1: Decimate by 50 using a 1 kHz cutoff FIR low-pass filter.
[iq, Fs] = etapa_decimacion(iq,Fs,1000,512,50); 

% Stage 2: Decimate by 20 using a 100 Hz cutoff FIR low-pass filter.
[iq, Fs] = etapa_decimacion(iq,Fs, 100,256,20);

% Stage 3: Decimate by 10 using a 20 Hz cutoff FIR low-pass filter.
[iq, Fs] = etapa_decimacion(iq,Fs,  20,256,10);

if abs(Fs-100) > 1e-9
    warning('Frontend output sampling rate is %.3f Hz, not 100 Hz.',Fs);
end

if debug_plots
    plot_spectrum(iq-mean(iq),Fs,8192,'Output FFT (DC removed just for illustration)','Hz',1);
end

fprintf('Frontend: %d samples output at %.1f Hz (%.1f s)\n', ...
    numel(iq),Fs,numel(iq)/Fs);
iq_out = iq(:);
end

% -------------------------------------------------------------------------
% ETAPA_DECIMACION Applies an FIR low-pass filter and downsamples
% the input signal by the specified decimation factor.
%
% Inputs:
%   x      - Input signal.
%   Fs_in  - Input sampling frequency (Hz).
%   fc     - FIR low-pass filter cutoff frequency (Hz).
%   order  - FIR filter order.
%   R      - Decimation factor.
%
% Outputs:
%   y      - Filtered and decimated signal.
%   Fs_out - Output sampling frequency (Hz).
% -------------------------------------------------------------------------
function [y, Fs_out] = etapa_decimacion(x, Fs_in, fc, order, R)
Fs_out = Fs_in/R;
nyq_out = Fs_out/2;

if fc >= Fs_in/2
    error('fc %.3f Hz must be less than Nyquist %.3f Hz.',fc,Fs_in/2);
end
if numel(x) <= 3*order
    error('Signal too short for FIR of order %d.',order);
end

transicion_aprox = 4*Fs_in/order;
if fc + transicion_aprox > nyq_out
    warning(['R=%d: The estimated FIR transition reaches the output Nyquist frequency. ' ...
        'It is a good idea to check for aliasing using freqz if you change Fs, order, or fc.'],R);
end

b = fir1(order,fc/(Fs_in/2));
yr = upfirdn(real(x),b,1,R);
yi = upfirdn(imag(x),b,1,R);
y = complex(yr,yi);

% Compensate for the FIR group delay and retain the expected
% number of output samples.
d = floor((order/2)/R);
n_out = floor(numel(x)/R);
i0 = d+1;
i1 = min(numel(y),i0+n_out-1);
y = y(i0:i1);
end

% -------------------------------------------------------------------------
% ESTIMAR_TONO Estimates the received Low-IF tone frequency.
%
% Uses a centered signal segment, Hamming window, and FFT-based
% peak search within a 5 kHz band around the nominal frequency.
% Parabolic interpolation is applied to refine the peak location.
% The spectral contrast is calculated relative to the local
% median magnitude within the search band.
%
% Outputs:
%   f_tone     - Estimated tone frequency (Hz).
%   contrast_db - Peak-to-local-median spectral contrast (dB).
% -------------------------------------------------------------------------
function [f_tone, contrast_db] = estimar_tono(x, Fs, f_nominal)
Nmax = min(numel(x),2^20);
i0 = max(1,floor((numel(x)-Nmax)/2)+1);
seg = x(i0:i0+Nmax-1);

Nfft = 2^nextpow2(Nmax);
X = fft(seg.*hamming(Nmax),Nfft);
f = (0:Nfft-1).'*Fs/Nfft;
f(f >= Fs/2) = f(f >= Fs/2)-Fs;

banda = abs(f-f_nominal) <= 5e3;
if ~any(banda)
    f_tone = f_nominal;
    contrast_db = -Inf;
    return;
end

% Calculate the signal-to-noise contrast (in dB) by comparing the peak amplitude to the local median
P = abs(X);
Psearch = P;
Psearch(~banda) = -Inf;
[~,k] = max(Psearch);

local_med = median(P(banda));
contrast_db = 20*log10((P(k)+eps)/(local_med+eps));

delta = 0;
if k > 1 && k < Nfft
    al = P(k-1); be = P(k); ga = P(k+1);
    den = al - 2*be + ga;
    if abs(den) > 1e-12
        delta = 0.5*(al-ga)/den;
        delta = max(-0.5,min(0.5,delta));
    end
end
f_tone = f(k) + delta*(Fs/Nfft);

offset_hz = f_tone - f_nominal;

if abs(offset_hz) > 2e3
    warning(['Low-IF deviation of %+.2f Hz from the nominal value. ' ...
        'The measured frequency for the digital mix is retained.'], ...
        offset_hz);
end
end

% -------------------------------------------------------------------------
% PLOT_SPECTRUM Displays the zero-centered FFT magnitude spectrum
% of the input signal for diagnostic purposes.
% -------------------------------------------------------------------------
function plot_spectrum(x, Fs, nfft_max, ttl, unit_label, unit_scale)
Nfft = min(nfft_max,numel(x));
X = fftshift(fft(x(1:Nfft)));
f = (-Nfft/2:Nfft/2-1)*(Fs/Nfft);
figure;
plot(f/unit_scale,20*log10(abs(X)+1e-12));
grid on; xlabel(unit_label); ylabel('dB'); title(ttl);
end
