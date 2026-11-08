%% ECHO - Estudio de ablacion (muestras estandar S1-S8)
% Ejecutar en la carpeta ECHO que contiene S1_1.mat ... S8_10.mat,
% frontend_dsp.m, procesar_signos_vitales_original.m y el archivo de
% referencias "Muestras ECHO(2).xlsx" (si está disponible).
%
% No modifica las capturas ni los procesadores original/experimental.
% La función procesar_signos_vitales_ablacion.m es una copia independiente
% del procesador original con tres interruptores de ablación.

clear; clc;

archivo_referencias = '';
candidatos_referencias = {'Muestras ECHO(2).xlsx','Muestras ECHO.xlsx','Muestras ECHO(1).xlsx'};
for ir = 1:numel(candidatos_referencias)
    if isfile(candidatos_referencias{ir})
        archivo_referencias = candidatos_referencias{ir};
        break;
    end
end
hoja_referencias = 'S 50cm-Ant 15cm';

%% 1. Detectar todas las capturas estándar existentes S1-S8
lista = dir('S*.mat');
archivos = {};
sujetos = [];
muestras = [];

for k = 1:numel(lista)
    tok = regexp(lista(k).name, '^S([1-8])_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tok), continue; end

    s = str2double(tok{1});
    m = str2double(tok{2});

    % Se incluyen todos los sufijos existentes del conjunto estándar,
    % incluidos _9 y _10 para S6-S8.
    if m >= 1 && m <= 10
        archivos{end+1,1} = lista(k).name; %#ok<SAGROW>
        sujetos(end+1,1) = s; %#ok<SAGROW>
        muestras(end+1,1) = m; %#ok<SAGROW>
    end
end

if isempty(archivos)
    error('No se encontraron capturas S1-S8 en la carpeta actual.');
end

[~, ord] = sortrows([sujetos, muestras]);
archivos = archivos(ord);
sujetos = sujetos(ord);
muestras = muestras(ord);
N = numel(archivos);

fprintf('Capturas estándar encontradas: %d\n', N);
if N ~= 70
    warning('Se esperaban 70 capturas estándar; se encontraron %d. Revise la carpeta.', N);
end

%% 2. Definir las variantes
nombres_variantes = { ...
    'Original_control', ...
    'Base_ablacion', ...
    'Sin_correccion_portadora', ...
    'Centro_IQ_estatico_Kasa', ...
    'Sin_detrend_cubico'};

nV = numel(nombres_variantes);

opciones = cell(nV,1);
for v = 1:nV
    opciones{v} = struct( ...
        'use_carrier_correction', true, ...
        'use_dynamic_center', true, ...
        'use_cubic_detrend', true);
end
% El primer control llama al procesador original, sin opciones.
% La segunda variante es la copia de ablación con todas las etapas activas.
opciones{3}.use_carrier_correction = false;
opciones{4}.use_dynamic_center = false;
opciones{5}.use_cubic_detrend = false;

%% 3. Cargar referencias reales si el Excel está disponible
RR_real_captura = nan(N,1);
HR_real_captura = nan(N,1);
if ~isempty(archivo_referencias) && isfile(archivo_referencias)
    [RR_real_captura, HR_real_captura] = cargar_referencias( ...
        archivo_referencias, hoja_referencias, sujetos, muestras);
    fprintf('Referencias cargadas desde: %s | hoja: %s\n', ...
        archivo_referencias, hoja_referencias);
else
    warning(['No se encontró "%s". Se ejecutará la ablación, pero las métricas ' ...
        'de error requerirán añadir referencias después.'], archivo_referencias);
end

%% 4. Reservar resultados en formato largo (una fila por captura y variante)
M = N*nV;
Archivo = cell(M,1);
Participante = nan(M,1);
Muestra = nan(M,1);
Variante = cell(M,1);
RR_real = nan(M,1);
HR_real = nan(M,1);
RR_rpm = nan(M,1);
HR_bpm = nan(M,1);
HR_candidato_bpm = nan(M,1);
HR_acuerdo = nan(M,1);
HR_baja_confianza = nan(M,1);
RR_baja_confianza = nan(M,1);
HR_prominencia = nan(M,1);
HR_score_ratio = nan(M,1);
HR_autocorr_bpm = nan(M,1);
HR_autocorr_pico = nan(M,1);
HR_espectral_bruto_bpm = nan(M,1);
HR_espectral_fundamental_bpm = nan(M,1);
HR_subarmonico_ratio = nan(M,1);
HR_consenso_bpm = nan(M,1);
HR_consenso_acuerdo = nan(M,1);
HR_consenso_dispersion = nan(M,1);
HR_fusion_count = nan(M,1);
RR_RMS_fase = nan(M,1);
HR_RMS_fase = nan(M,1);
Centro_offset_ratio = nan(M,1);
Arco_span_rad = nan(M,1);
Radio_circulo = nan(M,1);
Offset_portadora_Hz = nan(M,1);
Estado = repmat({''},M,1);
Error = repmat({''},M,1);

%% 5. Procesar cada captura con las mismas señales de frontend
for k = 1:N
    fprintf('\n[%d/%d] %s\n', k, N, archivos{k});

    try
        S = load(archivos{k});
        if isfield(S,'rx_continuous')
            rx = S.rx_continuous;
        elseif isfield(S,'rx')
            rx = S.rx;
        else
            error('El MAT no contiene rx_continuous ni rx.');
        end
        if ~isfield(S,'cfg')
            error('El MAT no contiene cfg.');
        end
        cfg = S.cfg;
        cfg.debug_plots = false;
        clear S

        % Frontend idéntico para todas las variantes de esta captura.
        [iq, Fs] = frontend_dsp(rx, cfg);
        clear rx

        dHP = designfilt('highpassiir', ...
            'FilterOrder',4, ...
            'HalfPowerFrequency',0.04, ...
            'SampleRate',Fs);
        dLP_rr = designfilt('lowpassiir', ...
            'FilterOrder',4, ...
            'HalfPowerFrequency',0.8, ...
            'SampleRate',Fs);

        for v = 1:nV
            r = (k-1)*nV + v;
            Archivo{r} = archivos{k};
            Participante(r) = sujetos(k);
            Muestra(r) = muestras(k);
            Variante{r} = nombres_variantes{v};
            RR_real(r) = RR_real_captura(k);
            HR_real(r) = HR_real_captura(k);

            try
                % El control original ejecuta la función original intacta.
                if v == 1
                    % evalc suprime la salida de diagnóstico extensa.
                    evalc('[rr, hr, info] = procesar_signos_vitales_original(iq, Fs, dHP, dLP_rr);');
                else
                    evalc(sprintf(['[rr, hr, info] = procesar_signos_vitales_ablacion(' ...
                        'iq, Fs, dHP, dLP_rr, opciones{%d});'], v));

                    % Para la variante base, las tres opciones están activas.
                end

                RR_rpm(r) = rr;
                HR_bpm(r) = hr;
                HR_candidato_bpm(r) = campo_num(info,'hr_candidate_bpm');
                HR_acuerdo(r) = campo_num(info,'hr_agreement');
                HR_baja_confianza(r) = campo_num(info,'hr_low_confidence');
                RR_baja_confianza(r) = campo_num(info,'rr_low_confidence');
                HR_prominencia(r) = campo_num(info,'hr_prominence');
                HR_score_ratio(r) = campo_num(info,'hr_score_ratio');
                HR_autocorr_bpm(r) = campo_num(info,'hr_autocorr_bpm');
                HR_autocorr_pico(r) = campo_num(info,'hr_autocorr_peak');
                HR_espectral_bruto_bpm(r) = campo_num(info,'hr_spectral_bpm_raw');
                HR_espectral_fundamental_bpm(r) = campo_num(info,'hr_spectral_bpm_fundamental');
                HR_subarmonico_ratio(r) = campo_num(info,'hr_subharmonic_ratio');
                HR_consenso_bpm(r) = campo_num(info,'hr_folded_bpm');
                HR_consenso_acuerdo(r) = campo_num(info,'hr_folded_agreement');
                HR_consenso_dispersion(r) = campo_num(info,'hr_folded_spread');
                HR_fusion_count(r) = campo_num(info,'hr_fusion_count');
                RR_RMS_fase(r) = campo_num(info,'rr_signal_rms_rad');
                HR_RMS_fase(r) = campo_num(info,'hr_signal_rms_rad');
                Centro_offset_ratio(r) = campo_num(info,'center_offset_ratio');
                Arco_span_rad(r) = campo_num(info,'arc_span_rad');
                Radio_circulo(r) = campo_num(info,'circle_radius');
                Offset_portadora_Hz(r) = campo_num(info,'carrier_offset_hz');

                if isfinite(hr), Estado{r} = 'HR_ACEPTADA';
                else, Estado{r} = 'HR_NaN'; end
            catch MEv
                Estado{r} = 'ERROR';
                Error{r} = MEv.message;
            end
        end

        clear iq dHP dLP_rr
    catch ME
        for v = 1:nV
            r = (k-1)*nV + v;
            Archivo{r} = archivos{k};
            Participante(r) = sujetos(k);
            Muestra(r) = muestras(k);
            Variante{r} = nombres_variantes{v};
            RR_real(r) = RR_real_captura(k);
            HR_real(r) = HR_real_captura(k);
            Estado{r} = 'ERROR_FRONTEND';
            Error{r} = ME.message;
        end
    end

    % Respaldo incremental para no perder todo si MATLAB se interrumpe.
    save('evaluar_ablacion_ECHO_progreso.mat', ...
        'archivos','sujetos','muestras','nombres_variantes', ...
        'Archivo','Participante','Muestra','Variante','RR_real','HR_real', ...
        'RR_rpm','HR_bpm','HR_candidato_bpm','HR_acuerdo', ...
        'HR_baja_confianza','RR_baja_confianza','HR_prominencia', ...
        'HR_score_ratio','HR_autocorr_bpm','HR_autocorr_pico', ...
        'HR_espectral_bruto_bpm','HR_espectral_fundamental_bpm', ...
        'HR_subarmonico_ratio','HR_consenso_bpm','HR_consenso_acuerdo', ...
        'HR_consenso_dispersion','HR_fusion_count','RR_RMS_fase','HR_RMS_fase', ...
        'Centro_offset_ratio','Arco_span_rad','Radio_circulo', ...
        'Offset_portadora_Hz','Estado','Error','k');
end

%% 6. Exportar resultados detallados
Resultados = table( ...
    Archivo, Participante, Muestra, Variante, RR_real, HR_real, RR_rpm, HR_bpm, ...
    HR_candidato_bpm, HR_acuerdo, HR_baja_confianza, RR_baja_confianza, ...
    HR_prominencia, HR_score_ratio, HR_autocorr_bpm, HR_autocorr_pico, ...
    HR_espectral_bruto_bpm, HR_espectral_fundamental_bpm, HR_subarmonico_ratio, ...
    HR_consenso_bpm, HR_consenso_acuerdo, HR_consenso_dispersion, HR_fusion_count, ...
    RR_RMS_fase, HR_RMS_fase, Centro_offset_ratio, Arco_span_rad, Radio_circulo, ...
    Offset_portadora_Hz, Estado, Error);

writetable(Resultados,'resultados_ablacion_ECHO.xlsx','Sheet','Resultados');
writetable(Resultados,'resultados_ablacion_ECHO.csv');

%% 7. Resumen de error y cobertura
VarianteResumen = cell(nV,1);
Capturas = zeros(nV,1);
HR_aceptadas = zeros(nV,1);
HR_rechazadas = zeros(nV,1);
HR_referencias = zeros(nV,1);
HR_cobertura_ref_pct = nan(nV,1);
HR_MAE_bpm = nan(nV,1);
HR_RMSE_bpm = nan(nV,1);
RR_aceptadas = zeros(nV,1);
RR_referencias = zeros(nV,1);
RR_cobertura_ref_pct = nan(nV,1);
RR_MAE_rpm = nan(nV,1);
RR_RMSE_rpm = nan(nV,1);
Errores_ejecucion = zeros(nV,1);

for v = 1:nV
    idx = strcmp(Variante,nombres_variantes{v});
    VarianteResumen{v} = nombres_variantes{v};
    Capturas(v) = nnz(idx);
    HR_aceptadas(v) = nnz(idx & isfinite(HR_bpm));
    HR_rechazadas(v) = nnz(idx & ~isfinite(HR_bpm) & ...
        ~strcmp(Estado,'ERROR') & ~strcmp(Estado,'ERROR_FRONTEND'));
    Errores_ejecucion(v) = nnz(idx & strcmp(Estado,'ERROR')) + nnz(idx & strcmp(Estado,'ERROR_FRONTEND'));

    hRef = idx & isfinite(HR_real);
    hEval = hRef & isfinite(HR_bpm);
    HR_referencias(v) = nnz(hRef);
    if nnz(hRef) > 0
        HR_cobertura_ref_pct(v) = 100*nnz(hEval)/nnz(hRef);
    end
    if any(hEval)
        e = HR_bpm(hEval)-HR_real(hEval);
        HR_MAE_bpm(v) = mean(abs(e));
        HR_RMSE_bpm(v) = sqrt(mean(e.^2));
    end

    rRef = idx & isfinite(RR_real);
    rEval = rRef & isfinite(RR_rpm);
    RR_referencias(v) = nnz(rRef);
    RR_aceptadas(v) = nnz(idx & isfinite(RR_rpm));
    if nnz(rRef) > 0
        RR_cobertura_ref_pct(v) = 100*nnz(rEval)/nnz(rRef);
    end
    if any(rEval)
        e = RR_rpm(rEval)-RR_real(rEval);
        RR_MAE_rpm(v) = mean(abs(e));
        RR_RMSE_rpm(v) = sqrt(mean(e.^2));
    end
end

Resumen = table(VarianteResumen, Capturas, HR_aceptadas, HR_rechazadas, ...
    HR_referencias, HR_cobertura_ref_pct, HR_MAE_bpm, HR_RMSE_bpm, ...
    RR_aceptadas, RR_referencias, RR_cobertura_ref_pct, RR_MAE_rpm, ...
    RR_RMSE_rpm, Errores_ejecucion);

writetable(Resumen,'resultados_ablacion_ECHO.xlsx','Sheet','Resumen');

fprintf('\n========== RESUMEN ABLACION ECHO ==========\n');
disp(Resumen);
fprintf('Archivos generados:\n');
fprintf('  resultados_ablacion_ECHO.xlsx\n');
fprintf('  resultados_ablacion_ECHO.csv\n');
fprintf('  evaluar_ablacion_ECHO_progreso.mat\n');

%% Funciones auxiliares locales
function valor = campo_num(info,campo)
    valor = NaN;
    if isstruct(info) && isfield(info,campo)
        x = info.(campo);
        if (isnumeric(x) || islogical(x)) && isscalar(x)
            valor = double(x);
        end
    end
end

function [rrReal, hrReal] = cargar_referencias(nombreArchivo, nombreHoja, sujetos, muestras)
    n = numel(sujetos);
    rrReal = nan(n,1);
    hrReal = nan(n,1);

    try
        C = readcell(nombreArchivo,'Sheet',nombreHoja);
    catch ME
        warning('No se pudieron leer referencias: %s', ME.message);
        return;
    end

    mapa = containers.Map('KeyType','char','ValueType','any');
    for i = 2:size(C,1)
        if size(C,2) < 3 || isempty(C{i,1}), continue; end
        id = strtrim(char(string(C{i,1})));
        tok = regexp(id,'^([1-8])_(\d+)$','tokens','once');
        if isempty(tok), continue; end

        rr = a_numero(C{i,2});
        hr = a_numero(C{i,3});
        mapa(id) = [rr,hr];
    end

    for k = 1:n
        id = sprintf('%d_%d',sujetos(k),muestras(k));
        if isKey(mapa,id)
            vals = mapa(id);
            rrReal(k) = vals(1);
            hrReal(k) = vals(2);
        end
    end
end

function x = a_numero(v)
    x = NaN;
    if isnumeric(v) && isscalar(v) && isfinite(v)
        x = double(v);
    elseif ischar(v) || isstring(v)
        y = str2double(strtrim(string(v)));
        if isscalar(y) && isfinite(y), x = double(y); end
    end
end
