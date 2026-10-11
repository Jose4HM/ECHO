function [iq_out, Fs] = frontend_dsp(iq_in, cfg)
%FRONTEND_DSP Low-IF a banda base para radar de signos vitales.
%
%   Entrada : IQ complejo crudo desde cw.m (Fs = cfg.fs, tono en cfg.tx_tone_hz)
%   Salida  : IQ complejo en banda base a 100 Hz, listo para
%             procesar_signos_vitales_fixed.
%
% -------------------------------------------------------------------------
% CAMBIOS v4:
%  [1] SALIDA A 100 Hz, no a 4 kHz. Las senales fisiologicas viven por
%      debajo de 3.5 Hz. A 4 kHz, un paso-alto de 0.1 Hz queda con polos a
%      radio 0.99992: mal condicionado, con transitorios de decenas de
%      segundos y ruido de redondeo justo en la banda cardiaca.
%      A 100 Hz todos los filtros aguas abajo quedan sanos.
%  [2] Se elimino `decimate`. Con r=50 en una sola etapa el Chebyshev
%      interno de orden 8 queda mal condicionado, MATLAB recomienda r<=13,
%      y `decimate` no esta pensado para entrada compleja. Ahora se usa
%      FIR de fase cero explicito + downsample en 3 etapas (50, 20, 10),
%      con la banda de rechazo garantizada por debajo del nuevo Nyquist en
%      cada paso.
%  [3] SE QUITO la resta de media DESPUES de mezclar. Ese termino DC de
%      banda base es la fuga TX->RX + clutter estatico, o sea el centro
%      del circulo IQ. Es informacion, no basura: el demodulador de fase
%      lo necesita. Solo se resta la media ANTES de mezclar, donde si
%      corresponde a DC real del ADC/LO a 0 Hz.
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% (a) Quitar el DC REAL del ADC / fuga de LO a 0 Hz.
%     Es legitimo: la senal util esta en +tx_tone_hz, asi que la media
%     sobre decenas de segundos solo captura el DC de hardware.
% -------------------------------------------------------------------------
iq = iq - mean(iq);

if debug_plots
    plot_spectrum(iq, Fs, 65536, 'Frontend v4 - FFT antes de mezclar', 'kHz', 1e3);
end

% -------------------------------------------------------------------------
% (b) Estimar la frecuencia REAL del tono low-IF.
%     TX y RX usan sintetizadores independientes (aunque compartan el
%     TCXO), asi que el tono recibido no cae exactamente en tx_tone_hz.
%     Mezclar con el valor nominal deja un residuo de frecuencia que
%     luego aparece como deriva rara en la fase.
% -------------------------------------------------------------------------
f_tone = estimar_tono(iq, Fs, cfg.tx_tone_hz);
fprintf('Frontend: tono nominal %.1f Hz, medido %.1f Hz (offset %.2f Hz)\n', ...
    cfg.tx_tone_hz, f_tone, f_tone - cfg.tx_tone_hz);

% -------------------------------------------------------------------------
% (c) Mezcla low-IF -> banda base con la frecuencia MEDIDA.
% -------------------------------------------------------------------------
n = (0:numel(iq)-1).';
iq = iq .* exp(-1j * 2*pi * f_tone * n / Fs);

% -------------------------------------------------------------------------
% (d) Decimacion en 3 etapas: 1 MHz -> 20 kHz -> 1 kHz -> 100 Hz.
%     En cada etapa la banda de rechazo del FIR queda por DEBAJO del nuevo
%     Nyquist, asi que no hay aliasing.
% -------------------------------------------------------------------------
[iq, Fs] = etapa_decimacion(iq, Fs, 1000, 512, 50);   % -> 20 kHz
[iq, Fs] = etapa_decimacion(iq, Fs,  100, 256, 20);   % -> 1 kHz
[iq, Fs] = etapa_decimacion(iq, Fs,   20, 256, 10);   % -> 100 Hz

% -------------------------------------------------------------------------
% (e) NO se resta la media aqui. Ver nota [3] arriba.
% -------------------------------------------------------------------------

if debug_plots
    plot_spectrum(iq - mean(iq), Fs, 8192, 'Frontend v4 - FFT salida (sin DC, solo para ver)', 'Hz', 1);
end

fprintf('Frontend: salida %d muestras a %.1f Hz (%.1f s)\n', ...
    numel(iq), Fs, numel(iq)/Fs);

iq_out = iq(:);
end

% =========================================================================
function [y, Fs_out] = etapa_decimacion(x, Fs_in, fc, order, R)
%ETAPA_DECIMACION FIR paso-bajo de fase cero + downsample por R.
%   Se verifica explicitamente que la banda de rechazo del FIR (aprox.
%   fc + 4*Fs_in/order) caiga por debajo del Nyquist de salida.

Fs_out = Fs_in / R;
nyq_out = Fs_out / 2;

transicion = 4 * Fs_in / order;   % ancho aproximado de la transicion FIR
f_stop = fc + transicion;
if f_stop > nyq_out
    warning(['Etapa de decimacion R=%d: la banda de rechazo empieza en ' ...
        '~%.1f Hz pero el Nyquist de salida es %.1f Hz. Sube el orden del ' ...
        'FIR o baja fc para evitar aliasing.'], R, f_stop, nyq_out);
end

if fc >= Fs_in/2
    error('fc=%.3f Hz debe ser menor que Nyquist=%.3f Hz.', fc, Fs_in/2);
end
if numel(x) <= 3*order
    error(['La senal tiene %d muestras, insuficientes para un FIR de ' ...
        'orden %d.'], numel(x), order);
end

bcoef = fir1(order, fc / (Fs_in/2));

% upfirdn filtra y decima en una sola pasada, calculando SOLO las muestras
% que sobreviven al downsampling: R veces mas rapido que filtrar todo y
% tirar el 98% de las muestras. Se aplica por separado a I y Q para no
% depender del soporte de entrada compleja.
yr = upfirdn(real(x), bcoef, 1, R);
yi = upfirdn(imag(x), bcoef, 1, R);
y = complex(yr, yi);

% Compensar el retardo de grupo del FIR (lineal, = order/2 muestras a la
% entrada) y recortar la cola de transitorio, para que el eje de tiempo de
% la salida siga correspondiendo al de la entrada.
d = floor((order/2) / R);
n_out = floor(numel(x) / R);
i0 = d + 1;
i1 = min(numel(y), i0 + n_out - 1);
y = y(i0:i1);
end

% =========================================================================
function f_tone = estimar_tono(x, Fs, f_nominal)
%ESTIMAR_TONO Localiza el pico espectral cerca de f_nominal e interpola.
%   Usa un tramo central de la captura para no depender de transitorios.

Nmax = min(numel(x), 2^20);
i0 = max(1, floor((numel(x) - Nmax)/2) + 1);
seg = x(i0:i0+Nmax-1);

Nfft = 2^nextpow2(Nmax);
X = fft(seg .* hamming(Nmax), Nfft);
f_ax = (0:Nfft-1).' * Fs / Nfft;
% Eje bilateral en [-Fs/2, Fs/2)
f_ax(f_ax >= Fs/2) = f_ax(f_ax >= Fs/2) - Fs;

% Buscar solo en un entorno de +-5 kHz alrededor del nominal
banda = abs(f_ax - f_nominal) <= 5e3;
if ~any(banda)
    warning('No se pudo acotar la busqueda del tono; se usa el valor nominal.');
    f_tone = f_nominal;
    return;
end

P = abs(X);
P(~banda) = -Inf;
[~, k] = max(P);

% Interpolacion parabolica para bajar del tamano de bin
if k > 1 && k < Nfft
    al = abs(X(k-1)); be = abs(X(k)); ga = abs(X(k+1));
    den = al - 2*be + ga;
    d = 0;
    if abs(den) > 1e-12
        d = 0.5 * (al - ga) / den;
        d = max(-0.5, min(0.5, d));
    end
else
    d = 0;
end

bin_hz = Fs / Nfft;
f_tone = f_ax(k) + d * bin_hz;

if abs(f_tone - f_nominal) > 2e3
    warning(['El tono medido (%.1f Hz) se aleja mas de 2 kHz del nominal ' ...
        '(%.1f Hz). Se usa el nominal por seguridad.'], f_tone, f_nominal);
    f_tone = f_nominal;
end
end

% =========================================================================
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
