% =========================================================================
% MASTER SCRIPT: RADAR DE SIGNOS VITALES
% =========================================================================
clear; clc; close all;

% --- PASO 1: CAPTURA (Tu código de bladeRF que ya funciona) ---

try
    [rx_continuous,cfg] = cw();
catch ME
    error(['Fallo en cw.m durante la captura: %s\n' ...
        'Revisa conexión del hardware y configuración antes de reintentar.'], ME.message);
end

if ~exist('rx_continuous','var')
    error('cw.m no generó rx_continuous.');
end

fprintf('\nEjecutando Front-End DSP...\n');

[iq_frontend,Fs_frontend] = frontend_dsp(rx_continuous,cfg);

% --- PASO 2: EL PUENTE (Decimación y preparación) ---
fs_raw = cfg.fs;      % 1e6 Hz
fs_proc = 100;        % 100 Hz (frecuencia ideal para signos vitales)
dec_factor = fs_raw / fs_proc;

if mod(dec_factor, 1) ~= 0
    error('fs_raw/fs_proc = %.4f no es entero; ajusta fs_proc para que dec_factor sea entero.', dec_factor);
end

% decimate() recomienda factores <= 13 por llamada (orden de filtro alto
% en una sola pasada es numéricamente inestable). Para dec_factor=10000
% esto es crítico: se deciman en cascada con factores pequeños en vez de
% una sola llamada directa.
dec_stages = local_factor_decimation(dec_factor, 13);
fprintf('Decimando en cascada: %s (factor total %d)\n', ...
    mat2str(dec_stages), dec_factor);

rx_decimado = rx_baseband;
for kstage = 1:numel(dec_stages)
    rx_decimado = decimate(rx_decimado, dec_stages(kstage));
end

% --- PASO 3: DISEÑO DE FILTROS ---
% Diseñamos los filtros que requiere tu función procesar_signos_vitales_fixed
dHP = designfilt('highpassiir', 'FilterOrder', 2, ...
                 'HalfPowerFrequency', 0.1, 'SampleRate', Fs_frontend);

dLP_rr = designfilt('lowpassiir', 'FilterOrder', 2, ...
                    'HalfPowerFrequency', 0.5, 'SampleRate', Fs_frontend);

% --- PASO 4: EJECUCIÓN DEL PIPELINE ---
fprintf('Datos capturados. Procesando signos vitales...\n');

[rpm,bpm] = procesar_signos_vitales_fixed_v2(iq_frontend, Fs_frontend, dHP, dLP_rr);

% --- PASO 5: RESULTADOS ---
fprintf('=============================\n');
fprintf('RESULTADOS FINALES:\n');
fprintf('Respiración: %.1f rpm\n', rpm);
fprintf('Ritmo Cardíaco: %.1f bpm\n', bpm);
fprintf('=============================\n');