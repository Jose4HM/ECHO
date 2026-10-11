function res = calibrar_posicion(seg)
%CALIBRAR_POSICION Verifica el montaje ANTES de medir, con el paciente sentado.
%
%   Captura corta (8 s por defecto), evalua los indicadores de geometria y
%   dice EXACTAMENTE que corregir. Repite hasta que salga APTO; solo
%   entonces hagas la medida real.
%
%   Los tres indicadores y su mando correspondiente:
%
%     Pico |rx|   -> 0.25-0.40   se ajusta con cfg.tx_amplitude
%     |c|/R       -> 2-10        se ajusta con la SEPARACION entre antenas
%     Excursion   -> 2-4 rad     depende del apuntado (altura y convergencia)
%
%   Son mandos casi independientes: la potencia mueve el nivel sin cambiar
%   la proporcion entre eco y fuga; la separacion cambia esa proporcion; el
%   apuntado cambia cuanta senal de pecho entra.
%
%   Uso:
%       calibrar_posicion        % 8 s
%       calibrar_posicion(12)    % 12 s, si la respiracion es muy lenta
clc;
if nargin < 1 || isempty(seg), seg = 8; end

fprintf('\n=== CALIBRACION DE POSICION ===\n');
fprintf('Paciente sentado, espalda apoyada, quieto, respirando normal.\n');
fprintf('Capturando %d s...\n\n', seg);

evalin('base', sprintf('cfg_cal_seg = %d;', seg));
[rx, cfg] = cw();
n = min(numel(rx), round(seg*cfg.fs));
rx = rx(1:n);

[iq, Fs] = frontend_dsp(rx, cfg);

dHP = designfilt('highpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',0.08,'SampleRate',Fs);
dLP = designfilt('lowpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',0.8,'SampleRate',Fs);

pico = max(abs(rx));
pct  = 100*nnz(abs(real(rx))>=0.99 | abs(imag(rx))>=0.99)/numel(rx);

% Geometria: ajuste de circulo global sobre el tramo corto
z = iq(:);
nn = (0:numel(z)-1)';
psi = unwrap(angle(z));
z = z .* exp(-1j*polyval(polyfit(nn, psi, 1), nn));   % quita rampa de LO
[c, R] = fit_circulo(z);
cR = abs(c)/max(R, eps);

Zc = z - c;
phi = unwrap(angle(Zc));
phi = phi - polyval(polyfit(nn, phi, 3), nn);
exc = prctile(filtfilt(dLP, phi), 97.5) - prctile(filtfilt(dLP, phi), 2.5);
arco = max(phi) - min(phi);
cerca = 100*nnz(abs(Zc) < 0.2*median(abs(Zc)))/numel(Zc);

res = struct('pico',pico,'pct_clip',pct,'cR',cR,'R',R, ...
             'excursion',exc,'arco',arco,'pct_origen',cerca);

fprintf('%-26s %8.3f   ', 'Pico |rx|', pico);
diagnostico(pico>=0.25 && pico<=0.40, pico<0.25, ...
    'sube cfg.tx_amplitude', 'baja cfg.tx_amplitude');

fprintf('%-26s %8.4f%%  ', 'Saturacion', pct);
diagnostico(pct < 0.001, false, '', 'baja cfg.tx_amplitude YA');

fprintf('%-26s %8.1f   ', '|c|/R', cR);
diagnostico(cR>=2 && cR<=10, cR<2, ...
    'ACERCA las antenas entre si', 'SEPARA las antenas entre si');

fprintf('%-26s %8.2f   ', 'Excursion respiratoria', exc);
diagnostico(exc>=2 && exc<=4, exc<2, ...
    'sube o baja la ALTURA hasta el 4o espacio intercostal', ...
    'el paciente respira muy fuerte: que respire normal');

fprintf('%-26s %8.2f   ', 'Apertura del arco (rad)', arco);
diagnostico(arco < 6, false, '', 'el paciente se mueve: insiste en inmovilidad');

fprintf('%-26s %8.1f%%  ', 'Muestras cerca del origen', cerca);
diagnostico(cerca < 2, false, '', 'ACERCA las antenas (|c|/R bajo)');

apto = pico>=0.25 && pico<=0.40 && pct<0.001 && cR>=1.5 && cR<=20 && ...
       exc>=1.5 && exc<=5 && arco<6 && cerca<3;

fprintf('\n');
if apto
    fprintf('>>> APTO. Puedes medir. NO muevas nada.\n');
else
    fprintf(2,'>>> NO APTO. Corrige lo marcado y repite calibrar_posicion.\n');
    fprintf(2,'    Orden recomendado: primero el nivel, luego |c|/R con la\n');
    fprintf(2,'    separacion, y la altura al final.\n');
end
res.apto = apto;
fprintf('===============================\n');
end

% =========================================================================
function diagnostico(ok, bajo, msg_bajo, msg_alto)
if ok
    fprintf('[OK]\n');
elseif bajo
    fprintf(2, '[X] %s\n', msg_bajo);
else
    fprintf(2, '[X] %s\n', msg_alto);
end
end

% =========================================================================
function [c, R] = fit_circulo(z)
x = real(z(:)); y = imag(z(:));
A = [x, y, ones(numel(x),1)];
b = x.^2 + y.^2;
p = A \ b;
c = complex(p(1)/2, p(2)/2);
R = sqrt(max(p(3) + abs(c)^2, eps));
end
