%% evaluar_compensacion_offset_ECHO.m
% Reproduce la evaluacion de compensacion de offset de frecuencia en ECHO.
%
% Entradas: archivos .mat con iq_frontend, Fs_frontend y vital_info.
%           vital_info.carrier_offset_hz debe contener el offset estimado
%           por el algoritmo original (en Hz).
% Salidas:  tabla de frecuencias equivalentes a pendientes de fase,
%           porcentaje de reduccion por captura, media y desviacion
%           estandar muestral. Opcionalmente exporta un CSV.
%
% Se reproduce el calculo usado para vacio_01.mat:
% pendiente antes  = -0.054173726 Hz
% pendiente despues = +0.001211391 Hz
% reduccion         = 97.7639 %%
%
% La prueba cuantifica reduccion de deriva LINEAL de fase, no reduccion
% del ruido de fase completo ni potencia RF recibida.

clearvars; clc;

%% CONFIGURACION
% Por defecto, se buscan los archivos junto a este script.
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta)
    carpeta = pwd;
end
patron = 'vacio_*.mat';

% Si prefieres nombres especificos, activa y ajusta esta lista:
archivos_explicitados = {};  % p. ej. {'vacio_01.mat','vacio_02.mat'}

inicio_s = 5;
fin_s = 55;
guardar_csv = true;
nombre_csv = 'resultados_offset_ECHO.csv';

%% IDENTIFICAR ARCHIVOS
if isempty(archivos_explicitados)
    listado = dir(fullfile(carpeta, patron));
    nombres = sort({listado.name});
else
    nombres = archivos_explicitados;
end

if isempty(nombres)
    error('No se encontraron archivos con el patron "%s" en %s', patron, carpeta);
end

N = numel(nombres);
archivo = strings(N,1);
offset_estimado_hz = nan(N,1);
pendiente_antes_hz = nan(N,1);
pendiente_despues_hz = nan(N,1);
reduccion_pct = nan(N,1);

fprintf('Evaluacion ECHO de compensacion de offset\n');
fprintf('Carpeta: %s\n', carpeta);
fprintf('Intervalo de regresion: %.1f a %.1f s\n', inicio_s, fin_s);
fprintf('Capturas encontradas: %d\n', N);

%% EVALUACION DE CADA CAPTURA
for k = 1:N
    archivo(k) = string(nombres{k});
    ruta = fullfile(carpeta, nombres{k});

    if ~isfile(ruta)
        error('No existe el archivo: %s', ruta);
    end

    % Solo se cargan variables pequenas. NO se carga rx_continuous.
    S = load(ruta, 'iq_frontend', 'Fs_frontend', 'vital_info');

    if ~isfield(S,'iq_frontend') || ~isfield(S,'Fs_frontend') || ...
       ~isfield(S,'vital_info') || ...
       ~isfield(S.vital_info,'carrier_offset_hz')
        error('Faltan variables requeridas en %s', nombres{k});
    end

    z = S.iq_frontend(:);
    Fs = double(S.Fs_frontend);
    df = double(S.vital_info.carrier_offset_hz);

    if isempty(z) || any(~isfinite(z)) || ~isfinite(Fs) || Fs <= 0 ...
            || ~isscalar(df) || ~isfinite(df)
        error('Datos invalidos en %s', nombres{k});
    end

    t = (0:numel(z)-1)' / Fs;
    idx = (t >= inicio_s) & (t <= fin_s);
    if nnz(idx) < 3
        error('No hay muestras suficientes entre %.1f y %.1f s en %s', ...
              inicio_s, fin_s, nombres{k});
    end

    % Se aplica la misma correccion registrada en el algoritmo ECHO.
    z_corregida = z .* exp(-1j * 2*pi * df * t);

    % Se ajusta una recta a la fase desenvolvida en el mismo intervalo.
    fase_antes = unwrap(angle(z));
    fase_despues = unwrap(angle(z_corregida));
    ajuste_antes = polyfit(t(idx), fase_antes(idx), 1);
    ajuste_despues = polyfit(t(idx), fase_despues(idx), 1);

    f_antes = ajuste_antes(1) / (2*pi);
    f_despues = ajuste_despues(1) / (2*pi);

    if abs(f_antes) <= eps(max(1, abs(df)))
        error('Pendiente previa casi nula en %s: reduccion indefinida.', nombres{k});
    end

    porcentaje = 100 * (1 - abs(f_despues)/abs(f_antes));

    offset_estimado_hz(k) = df;
    pendiente_antes_hz(k) = f_antes;
    pendiente_despues_hz(k) = f_despues;
    reduccion_pct(k) = porcentaje;

    fprintf('\n===== %s =====\n', nombres{k});
    fprintf('Offset estimado:   %+.9f Hz\n', df);
    fprintf('Pendiente antes:   %+.9f Hz\n', f_antes);
    fprintf('Pendiente despues: %+.9f Hz\n', f_despues);
    fprintf('Reduccion:        %.4f %%\n', porcentaje);
end

%% RESUMEN
resultados = table(archivo, offset_estimado_hz, pendiente_antes_hz, ...
    pendiente_despues_hz, reduccion_pct, ...
    'VariableNames', {'Captura','OffsetEstimado_Hz','PendienteAntes_Hz', ...
                      'PendienteDespues_Hz','Reduccion_pct'});

fprintf('\n===== TABLA DE RESULTADOS =====\n');
disp(resultados);

fprintf('Media de reduccion: %.4f %%\n', mean(reduccion_pct));
if N >= 2
    fprintf('Desviacion estandar muestral: %.4f puntos porcentuales\n', ...
        std(reduccion_pct, 0));
else
    fprintf('Desviacion estandar muestral: no definida (N = 1)\n');
end
fprintf('Minimo: %.4f %% | Maximo: %.4f %% | N = %d\n', ...
        min(reduccion_pct), max(reduccion_pct), N);

%% EXPORTACION OPCIONAL
if guardar_csv
    destino = fullfile(carpeta, nombre_csv);
    writetable(resultados, destino);
    fprintf('CSV guardado en: %s\n', destino);
end
