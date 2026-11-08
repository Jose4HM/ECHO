function [bpm_ac, detalle] = hr_autocorrelacion(info, hr_min, hr_max)
%HR_AUTOCORRELACION Estimador de ritmo cardiaco por autocorrelacion.
%
%   Alternativa al pico espectral, pensada para el caso en que el latido
%   compite con armonicos altos de la respiracion.
%
%   IDEA: en el espectro, un armonico en 5*RR y un latido a esa misma
%   frecuencia son la misma raya. En la autocorrelacion no: los armonicos
%   4x, 5x, 6x pertenecen a UNA sola senal periodica cuyo periodo es el de
%   la respiracion, y refuerzan la correlacion en ~4 s, no en 0.7 s. El
%   latido si es una periodicidad propia en 0.5-1.1 s.
%
%   Ademas se penalizan explicitamente los retardos que son submultiplos
%   del periodo respiratorio, que es donde caen los armonicos aislados.
%
%   Uso:
%       [~, ~, info] = replotear('captura.mat');
%       hr_autocorrelacion(info)

if nargin < 2 || isempty(hr_min), hr_min = 55;  end
if nargin < 3 || isempty(hr_max), hr_max = 115; end

x  = info.x_hr(:);
Fs = info.Fs;
x  = x - mean(x);

% Autocorrelacion normalizada, solo retardos positivos
r = xcorr(x, 'coeff');
r = r(numel(x):end);
lags_s = (0:numel(r)-1)'/Fs;

lag_lo = round(60/hr_max * Fs);
lag_hi = round(60/hr_min * Fs);
lag_hi = min(lag_hi, numel(r)-2);
idx = (lag_lo:lag_hi)';

score = r(idx);

% Penalizacion de submultiplos del periodo respiratorio. Un armonico
% aislado en h*f_rr produce un pico en T_rr/h; el latido no tiene por que
% caer ahi. Se atenuan esos retardos un 40%.
T_rr = 1 / max(info.f_rr_est, 1e-6);
pen  = ones(size(idx));
for h = 3:10
    T_h = T_rr / h;
    if T_h >= lags_s(lag_lo) && T_h <= lags_s(lag_hi)
        pen = pen .* (1 - 0.40 * exp(-((lags_s(idx) - T_h).^2) / (2*(0.020)^2)));
    end
end
score_pen = score .* pen;

[~, k] = max(score_pen);
kk = idx(k);

% Interpolacion parabolica sobre el retardo
if kk > 1 && kk < numel(r)-1
    y0 = r(kk-1); y1 = r(kk); y2 = r(kk+1);
    den = y0 - 2*y1 + y2;
    d = 0;
    if abs(den) > 1e-12, d = 0.5*(y0-y2)/den; d = max(-0.5,min(0.5,d)); end
else
    d = 0;
end
lag_ref = (kk-1) + d;
bpm_ac  = 60 * Fs / lag_ref;

detalle.r = r;
detalle.lags_s = lags_s;
detalle.pico_correlacion = r(kk);
detalle.T_rr = T_rr;

fprintf('\n--- HR POR AUTOCORRELACION ---\n');
fprintf('Periodo respiratorio  : %.2f s (%.1f rpm)\n', T_rr, 60/T_rr);
fprintf('Retardo del pico      : %.3f s\n', lag_ref/Fs);
fprintf('Correlacion en el pico: %.3f\n', r(kk));
fprintf('HR estimado           : %.1f bpm\n', bpm_ac);
if r(kk) < 0.25
    fprintf(2,'[!] Correlacion baja (<0.25): periodicidad debil, poco fiable.\n');
end
fprintf('------------------------------\n');

figure('Name','HR por autocorrelacion');
subplot(2,1,1);
plot(lags_s(1:min(end,round(8*Fs))), r(1:min(end,round(8*Fs))), 'LineWidth',1.1);
hold on; grid on;
xline(T_rr, '--', 'periodo respiratorio', 'Color',[0.6 0.6 0.6]);
xline(lag_ref/Fs, '--r', sprintf('%.1f bpm', bpm_ac));
xlabel('retardo (s)'); ylabel('correlacion');
title('Autocorrelacion de la banda cardiaca');

subplot(2,1,2);
plot(60./lags_s(idx), score, 'LineWidth',1.1); hold on;
plot(60./lags_s(idx), score_pen, '--', 'LineWidth',1.1);
grid on; xlabel('bpm'); ylabel('correlacion');
legend('sin penalizar','penalizando submultiplos de RR','Location','best');
xline(bpm_ac, '--r');
title('Puntuacion en la banda cardiaca');
end
