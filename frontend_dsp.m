function [iq_out, Fs] = frontend_dsp(iq_in, cfg)
%FRONTEND_DSP Convierte el tono Low-IF recibido a banda base y 100 Hz.

if nargin < 2 || ~isstruct(cfg) || ~isfield(cfg,'fs') || ~isfield(cfg,'tx_tone_hz')
    error('cfg debe contener fs y tx_tone_hz.');
end

Fs = cfg.fs;
if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('cfg.fs debe ser positivo.');
end

iq = iq_in(:);
if isempty(iq)
    error('iq_in esta vacio.');
end

if any(~isfinite(iq))
    good = isfinite(iq);
    if nnz(good) < 2
        error('iq_in no tiene suficientes muestras finitas.');
    end
    idx = (1:numel(iq)).';
    iq = interp1(idx(good),iq(good),idx,'linear','extrap');
end

debug_plots = isfield(cfg,'debug_plots') && logical(cfg.debug_plots);

% El DC a 0 Hz antes de mezclar es offset del receptor. La senal radar util
% esta alrededor de +tx_tone_hz, por lo que esta resta no elimina clutter
% de banda base; ese clutter aparece recien despues de la mezcla.
iq = iq - mean(iq);

if debug_plots
    plot_spectrum(iq,Fs,65536,'FFT antes de mezclar','kHz',1e3);
end

[f_tone, tone_contrast_db] = estimar_tono(iq,Fs,cfg.tx_tone_hz);
fprintf('Frontend: tono nominal %.1f Hz, medido %.2f Hz (offset %+.2f Hz)\n', ...
    cfg.tx_tone_hz,f_tone,f_tone-cfg.tx_tone_hz);
fprintf('Frontend: contraste del tono %.1f dB sobre la mediana local\n',tone_contrast_db);
if tone_contrast_db < 10
    warning(['El tono Low-IF tiene solo %.1f dB de contraste. Revise TX, ' ...
        'antenas, interferencia o saturacion antes de confiar en RR/HR.'], ...
        tone_contrast_db);
end

% Mezcla a banda base con la frecuencia realmente observada.
n = (0:numel(iq)-1).';
iq = iq .* exp(-1j*2*pi*f_tone*n/Fs);

% Decimacion multietapa 1 MHz -> 20 kHz -> 1 kHz -> 100 Hz.
[iq, Fs] = etapa_decimacion(iq,Fs,1000,512,50);
[iq, Fs] = etapa_decimacion(iq,Fs, 100,256,20);
[iq, Fs] = etapa_decimacion(iq,Fs,  20,256,10);

if abs(Fs-100) > 1e-9
    warning('La salida del frontend es %.3f Hz, no 100 Hz.',Fs);
end

% No quitar la media aqui: ahora el DC contiene fuga/clutter y define el
% centro de la trayectoria I/Q que usa el demodulador de fase.
if debug_plots
    plot_spectrum(iq-mean(iq),Fs,8192,'FFT de salida (DC quitado solo para ver)','Hz',1);
end

fprintf('Frontend: salida %d muestras a %.1f Hz (%.1f s)\n', ...
    numel(iq),Fs,numel(iq)/Fs);
iq_out = iq(:);
end

function [y, Fs_out] = etapa_decimacion(x, Fs_in, fc, order, R)
Fs_out = Fs_in/R;
nyq_out = Fs_out/2;

if fc >= Fs_in/2
    error('fc %.3f Hz debe ser menor que Nyquist %.3f Hz.',fc,Fs_in/2);
end
if numel(x) <= 3*order
    error('Senal demasiado corta para FIR de orden %d.',order);
end

% Comprobacion conservadora de la transicion aproximada de una Hamming.
transicion_aprox = 4*Fs_in/order;
if fc + transicion_aprox > nyq_out
    warning(['R=%d: la transicion FIR estimada alcanza el Nyquist de salida. ' ...
        'Conviene revisar aliasing con freqz si se cambian Fs/order/fc.'],R);
end

b = fir1(order,fc/(Fs_in/2));
yr = upfirdn(real(x),b,1,R);
yi = upfirdn(imag(x),b,1,R);
y = complex(yr,yi);

% Compensacion del retardo de grupo de la FIR lineal.
d = floor((order/2)/R);
n_out = floor(numel(x)/R);
i0 = d+1;
i1 = min(numel(y),i0+n_out-1);
y = y(i0:i1);
end

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

P = abs(X);
Psearch = P;
Psearch(~banda) = -Inf;
[~,k] = max(Psearch);

local_med = median(P(banda));
contrast_db = 20*log10((P(k)+eps)/(local_med+eps));

% Interpolacion parabolica alrededor del maximo.
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
if abs(f_tone-f_nominal) > 2e3
    warning('Tono medido %.1f Hz demasiado lejos del nominal; se usa %.1f Hz.', ...
        f_tone,f_nominal);
    f_tone = f_nominal;
end
end

function plot_spectrum(x, Fs, nfft_max, ttl, unit_label, unit_scale)
Nfft = min(nfft_max,numel(x));
X = fftshift(fft(x(1:Nfft)));
f = (-Nfft/2:Nfft/2-1)*(Fs/Nfft);
figure;
plot(f/unit_scale,20*log10(abs(X)+1e-12));
grid on; xlabel(unit_label); ylabel('dB'); title(ttl);
end
