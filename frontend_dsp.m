function [iq_out, Fs] = frontend_dsp(iq_in,cfg)
%% ==========================================================
% 1. ELIMINAR DC
% ===========================================================

iq = iq_in - mean(iq_in);

Fs = cfg.fs;

%% ==========================================================
% 2. FFT DE DIAGNÓSTICO
% ===========================================================

Nfft = min(65536,length(iq));

X = fftshift(fft(iq(1:Nfft)));

f = (-Nfft/2:Nfft/2-1)*(Fs/Nfft);

figure
plot(f/1000,20*log10(abs(X)+1e-12))
grid on

xlabel('kHz')
ylabel('dB')

title('Frontend - FFT antes de mezclar')

%% ==========================================================
% 3. MEZCLA A BANDA BASE
% ===========================================================

n = (0:length(iq)-1).';

iq = iq .* exp(-1j*2*pi*cfg.tx_tone_hz*n/Fs);

%% ==========================================================
% 4. LPF
% ===========================================================

fc = 20000;

N = 128;

b = fir1(N,fc/(Fs/2));

iq = filtfilt(b,1,iq);

%% ==========================================================
% 5. DECIMACIÓN (1MHz a 20Khz)
% ===========================================================

R1 = 50;

iq = decimate(iq,R1);

Fs = Fs/R1;

fc = 200;

N = 128;

%% ==========================================================
% 5.1 DECIMACIÓN (20Khz a 1000 Hz)
% ===========================================================

b = fir1(N,fc/(Fs/2));

iq = filtfilt(b,1,iq);

R2 = 20;

iq = decimate(iq,R2);

Fs = Fs/R2;

Nfft = min(8192,length(iq));

X = fftshift(fft(iq(1:Nfft)));

f = (-Nfft/2:Nfft/2-1)*(Fs/Nfft);

figure
plot(f,20*log10(abs(X)+1e-12))
grid on

xlabel('Hz')
ylabel('dB')

title('Frontend - FFT salida')

iq_out = iq;