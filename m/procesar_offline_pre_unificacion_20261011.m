function [rpm, bpm, info] = procesar_offline(archivo, rango_s)
%PROCESAR_OFFLINE Reprocesa una captura guardada, sin tocar el hardware.
%
%   Cada captura cuesta 60 s de bladeRF. Para iterar sobre el DSP no hace
%   falta volver a capturar: guarda una vez y reprocesa cuantas veces
%   quieras.
%
%   Para GUARDAR una captura, despues de correr main.m:
%       save('captura_01.mat', 'rx_continuous', 'cfg', '-v7.3');
%
%   Para REPROCESARLA:
%       [rpm, bpm, info] = procesar_offline('captura_01.mat');
%
%   Para procesar SOLO UN TRAMO (segundos [inicio fin]):
%       procesar_offline('apnea.mat', [0 17]);
%
%   El recorte es imprescindible en capturas de apnea: si la apnea dura
%   17 s y la captura 60, los 43 s de respiracion posterior dominan por
%   completo el espectro y anulan la ventaja de haber hecho la apnea.
%   Tambien sirve para saltarse un movimiento corporal puntual.
%
%   Recomendacion para la presentacion: guarda 8-10 capturas etiquetadas
%   con la verdad de referencia en el nombre, por ejemplo
%   'sujeto1_reposo_HR72_RR14.mat'. Con eso puedes hacer una diapositiva
%   de estimado vs. real sin depender de que el radar funcione en vivo.

if nargin < 1
    error('Uso: procesar_offline(''captura_01.mat'')');
end
if exist(archivo, 'file') ~= 2
    error('No se encontro el archivo %s', archivo);
end

S = load(archivo);
if ~isfield(S, 'rx_continuous') || ~isfield(S, 'cfg')
    error('El .mat debe contener las variables rx_continuous y cfg.');
end

rx_continuous = S.rx_continuous;
cfg = S.cfg;

dur_total = numel(rx_continuous)/cfg.fs;

if nargin >= 2 && ~isempty(rango_s)
    assert(numel(rango_s)==2 && rango_s(2) > rango_s(1), ...
        'rango_s debe ser [inicio fin] en segundos.');
    dur_tramo = min(rango_s(2), dur_total) - rango_s(1);

    if dur_tramo < 12
        warning(['Tramo de solo %.1f s. La resolucion espectral sera de ' ...
            '%.1f bpm y el recorte de bordes se comera buena parte. ' ...
            'Por debajo de ~12 s el numero no es fiable.'], ...
            dur_tramo, 60/dur_tramo);
    end

    % OJO: no se trunca rx_continuous. El rango se pasa al procesador, que
    % lo aplica DESPUES de estimar la deriva de LO y el centro del circulo
    % con la captura completa. Truncar aqui dejaria esas dos estimaciones
    % sin datos suficientes.
    fprintf('Reprocesando %s, analisis en %.1f-%.1f s de %.1f s totales (%.1f s)...\n', ...
        archivo, rango_s(1), rango_s(2), dur_total, dur_tramo);
else
    fprintf('Reprocesando %s (%.1f s a %.0f Hz)...\n', ...
        archivo, dur_total, cfg.fs);
end

[iq_frontend, Fs_frontend] = frontend_dsp(rx_continuous, cfg);

dHP = designfilt('highpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 0.08, 'SampleRate', Fs_frontend);
dLP_rr = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 0.8, 'SampleRate', Fs_frontend);

if nargin >= 2 && ~isempty(rango_s)
    [rpm, bpm, info] = procesar_signos_vitales_fixed( ...
        iq_frontend, Fs_frontend, dHP, dLP_rr, rango_s);
else
    [rpm, bpm, info] = procesar_signos_vitales_fixed( ...
        iq_frontend, Fs_frontend, dHP, dLP_rr);
end

fprintf('\n=============================\n');
if info.rr_fallback
    fprintf('Respiracion: SIN DETECCION\n');
else
    fprintf('Respiracion: %.1f rpm\n', rpm);
end
if info.hr_fallback
    fprintf('Ritmo cardiaco: SIN DETECCION\n');
else
    fprintf('Ritmo cardiaco: %.1f bpm (acuerdo %.0f%%)\n', bpm, 100*info.hr_agreement);
end
fprintf('=============================\n');
end
