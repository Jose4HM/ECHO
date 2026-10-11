function espectrograma_cardiaco(archivo, hr_real)
%ESPECTROGRAMA_CARDIACO Responde a UNA pregunta: ¿hay latido en los datos?
%
%   No estima nada ni decide nada. Dibuja la energia en la banda cardiaca
%   en funcion del tiempo y deja que lo veas.
%
%   Uso:
%       espectrograma_cardiaco('sujeto1_normal_HR80.mat', 80)
%
%   COMO LEERLO:
%
%   [A] Hay una BANDA HORIZONTAL continua cerca de la linea blanca.
%       -> El latido esta en los datos. El problema es el estimador y
%          merece la pena seguir con DSP.
%
%   [B] Hay trozos de banda que aparecen y desaparecen.
%       -> Senal marginal. El limite es fisico: geometria, apuntado o
%          distancia. Mas DSP no lo va a arreglar.
%
%   [C] Manchas dispersas sin estructura horizontal.
%       -> No hay latido detectable. Hay que cambiar el montaje antes de
%          seguir tocando codigo.
%
%   Se dibuja SIN cancelacion de respiracion, a proposito: queremos ver
%   los datos crudos, no el resultado de un algoritmo que ya sabemos que
%   puede inventarse cosas.

if nargin < 1, error('Uso: espectrograma_cardiaco(''captura.mat'', hr_real)'); end
if nargin < 2, hr_real = []; end

[~, ~, info] = procesar_offline(archivo);

x  = info.x(:);
Fs = 100;

% Banda cardiaca generosa: 42 a 150 bpm. Sin enmascarar armonicos.
d = designfilt('bandpassiir', 'FilterOrder', 6, ...
    'HalfPowerFrequency1', 0.70, 'HalfPowerFrequency2', 2.50, ...
    'SampleRate', Fs);
xc = filtfilt(d, x);

win_s   = 12;                     % resolucion 5 bpm
L       = round(win_s*Fs);
solape  = round(0.92*L);          % avance de ~1 s: traza suave
Nfft    = 2^15;

figure('Name', sprintf('Espectrograma cardiaco - %s', archivo));

subplot(2,1,1);
t = (0:numel(xc)-1)'/Fs;
plot(t, xc, 'LineWidth', 0.8); grid on;
xlabel('s'); ylabel('rad');
title('Banda cardiaca 0.7-2.5 Hz (sin cancelar respiracion)');

subplot(2,1,2);
[S, f, tt] = spectrogram(xc, hann(L), solape, Nfft, Fs);
bpm = f*60;
sel = bpm >= 40 & bpm <= 150;

P = abs(S(sel,:)).^2;
% Normalizar cada columna: interesa la FORMA del espectro en cada instante,
% no que un tramo tenga mas energia total que otro.
P = P ./ (max(P, [], 1) + eps);

imagesc(tt, bpm(sel), 10*log10(P + 1e-12));
axis xy; colormap(jet); caxis([-12 0]);
colorbar; xlabel('s'); ylabel('bpm');
hold on;
if ~isempty(hr_real)
    yline(hr_real, 'w--', sprintf('%.0f bpm real', hr_real), ...
        'LineWidth', 1.8, 'LabelHorizontalAlignment', 'left');
end
title('Energia normalizada por instante -- busca una BANDA HORIZONTAL');

fprintf('\n=== COMO LEER LA FIGURA ===\n');
fprintf('[A] Banda horizontal continua cerca de la linea blanca\n');
fprintf('    -> el latido esta ahi. Seguir con DSP tiene sentido.\n');
fprintf('[B] Trozos de banda intermitentes\n');
fprintf('    -> senal marginal. El limite es fisico, no de algoritmo.\n');
fprintf('[C] Manchas sin estructura\n');
fprintf('    -> no hay latido detectable. Cambiar montaje antes que codigo.\n');
end
