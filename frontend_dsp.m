function [iq_out, Fs] = frontend_dsp(iq_in, cfg)
%FRONTEND_DSP Low-IF to baseband frontend for vital-sign radar.
%   Input is raw complex IQ from cw.m. Output is complex baseband IQ at
%   4 kHz, ready for procesar_signos_vitales_v3.

if nargin < 2 || ~isstruct(cfg)
    error('cfg debe ser una estructura valida.');
end
if ~isfield(cfg, 'fs') || ~isfield(cfg, 'tx_tone_hz')
    error('cfg debe contener fs y tx_tone_hz.');
end

Fs = cfg.fs;
if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('cfg.fs debe ser un escalar positivo en Hz.');
end

iq = iq_in(:);
if isempty(iq)
    error('iq_in esta vacio.');
end

if any(~isfinite(iq))
    n_bad = sum(~isfinite(iq));
    warning('iq_in contiene %d muestras no finitas; se interpolaran.', n_bad);
    idx_good = isfinite(iq);
    if sum(idx_good) < 2
        error('iq_in no tiene suficientes muestras finitas.');
    end
    idx_all = (1:numel(iq)).';
    iq = interp1(idx_all(idx_good), iq(idx_good), idx_all, 'linear', 'extrap');
end

debug_plots = isfield(cfg, 'debug_plots') && logical(cfg.debug_plots);

iq = iq - mean(iq);

if debug_plots
    plot_spectrum(iq, Fs, 65536, 'Frontend v3 - FFT antes de mezclar', 'kHz', 1e3);
end

% Mezcla Low-IF a banda base.
n = (0:numel(iq)-1).';
iq = iq .* exp(-1j * 2*pi * cfg.tx_tone_hz * n / Fs);

% 1 MHz -> 100 kHz. Primera etapa de decimacion: deja margen amplio y
% evita pedir un filtro irrealmente angosto antes de bajar a 20 kHz.
iq = lowpass_fir_zero_phase(iq, Fs, 40000, 256);
iq = downsample(iq, 10);
Fs = Fs / 10;

% 100 kHz -> 20 kHz. Antes de decimar x5, Nyquist final sera 10 kHz;
% por eso el corte debe quedar por debajo, no en 20 kHz.
iq = lowpass_fir_zero_phase(iq, Fs, 8000, 256);
iq = downsample(iq, 5);
Fs = Fs / 5;

% 20 kHz -> 4 kHz. La banda de interes fisiologica queda muy por debajo
% de 100 Hz, pero se conserva margen para diagnostico y filtrado posterior.
iq = lowpass_fir_zero_phase(iq, Fs, 100, 256);
iq = downsample(iq, 5);
Fs = Fs / 5;

iq = iq - mean(iq);

if debug_plots
    plot_spectrum(iq, Fs, 8192, 'Frontend v3 - FFT salida', 'Hz', 1);
end

iq_out = iq(:);
end

function y = lowpass_fir_zero_phase(x, Fs, fc, order)
if fc >= Fs/2
    error('fc=%.3f Hz debe ser menor que Nyquist=%.3f Hz.', fc, Fs/2);
end
b = fir1(order, fc / (Fs/2));
y = filtfilt(b, 1, x);
end

function plot_spectrum(x, Fs, nfft_max, ttl, unit_label, unit_scale)
Nfft = min(nfft_max, numel(x));
X = fftshift(fft(x(1:Nfft)));
f = (-Nfft/2:Nfft/2-1) * (Fs / Nfft);
figure;
plot(f / unit_scale, 20*log10(abs(X) + 1e-12));
grid on;
xlabel(unit_label);
ylabel('dB');
title(ttl);
end
