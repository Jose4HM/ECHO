% =========================================================================
% MASTER SCRIPT: RADAR DE SIGNOS VITALES
% =========================================================================
clear; clc; close all;

% --- PASO 1: CAPTURA (Tu código de bladeRF que ya funciona) ---

run('cw.m');

% Alternativa segura si cw.m devuelve variables en su propio scope, usar feval de función
% (descomentar si cw.m está implementado como función cw)
% [rx_baseband, cfg] = feval('cw'); 
% --- PASO 2: EL PUENTE (Decimación y preparación) ---
fs_raw = cfg.fs;      % 1e6 Hz
fs_proc = 100;        % 100 Hz (frecuencia ideal para signos vitales)
dec_factor = fs_raw / fs_proc;

% Decimamos la señal: eliminamos altas frecuencias y reducimos tasa
rx_decimado = decimate(rx_baseband, dec_factor);

% --- PASO 3: DISEÑO DE FILTROS ---
% Diseñamos los filtros que requiere tu función procesar_signos_vitales_fixed
dHP = designfilt('highpassiir', 'FilterOrder', 2, ...
                 'HalfPowerFrequency', 0.1, 'SampleRate', fs_proc);

dLP_rr = designfilt('lowpassiir', 'FilterOrder', 2, ...
                    'HalfPowerFrequency', 0.5, 'SampleRate', fs_proc);

% --- PASO 4: EJECUCIÓN DEL PIPELINE ---
fprintf('Datos capturados. Procesando signos vitales...\n');

[rpm, bpm] = procesar_signos_vitales_fixed(rx_decimado, fs_proc, dHP, dLP_rr);

% --- PASO 5: RESULTADOS ---
fprintf('=============================\n');
fprintf('RESULTADOS FINALES:\n');
fprintf('Respiración: %.1f rpm\n', rpm);
fprintf('Ritmo Cardíaco: %.1f bpm\n', bpm);
fprintf('=============================\n');