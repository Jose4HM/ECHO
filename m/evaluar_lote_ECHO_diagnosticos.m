
%% ECHO - Comparacion de algoritmos (muestras estandar)
clear; clc;

% Solo participantes S1 a S8, capturas 1 a 10.
lista = dir('S*.mat');

archivos = {};
sujetos = [];
muestras = [];

for k = 1:numel(lista)
    t = regexp(lista(k).name, ...
        '^S([1-8])_(\d+)\.mat$', 'tokens', 'once');

    if isempty(t)
        continue;
    end

    s = str2double(t{1});
    m = str2double(t{2});

    if m >= 1 && m <= 10
        archivos{end+1,1} = lista(k).name;
        sujetos(end+1,1) = s;
        muestras(end+1,1) = m;
    end
end

% Ordenar por participante y numero de captura
[~, orden] = sortrows([sujetos, muestras]);
archivos = archivos(orden);
sujetos = sujetos(orden);
muestras = muestras(orden);

N = numel(archivos);
fprintf('Capturas estandar encontradas: %d\n', N);

assert(N > 0, 'No se encontraron capturas.');

if N ~= 70
    warning('Se esperaban 70 capturas. Se encontraron %d.', N);
end

% Resultados
RR_original = nan(N,1);
HR_original = nan(N,1);
RR_experimental = nan(N,1);
HR_experimental = nan(N,1);

Acuerdo_original = nan(N,1);
Acuerdo_experimental = nan(N,1);

Candidato_original = nan(N,1);
Candidato_experimental = nan(N,1);

Confianza_original = nan(N,1);
Confianza_experimental = nan(N,1);

% Diagnosticos internos de HR (se exportan aunque HR sea rechazada).
HR_prominencia_original = nan(N,1);
HR_prominencia_experimental = nan(N,1);
HR_score_ratio_original = nan(N,1);
HR_score_ratio_experimental = nan(N,1);
HR_espectral_bruto_original = nan(N,1);
HR_espectral_bruto_experimental = nan(N,1);
HR_espectral_fundamental_original = nan(N,1);
HR_espectral_fundamental_experimental = nan(N,1);
HR_ratio_subarmonico_original = nan(N,1);
HR_ratio_subarmonico_experimental = nan(N,1);
HR_autocorr_bpm_original = nan(N,1);
HR_autocorr_bpm_experimental = nan(N,1);
HR_autocorr_pico_original = nan(N,1);
HR_autocorr_pico_experimental = nan(N,1);
HR_consenso_ventanas_original = nan(N,1);
HR_consenso_ventanas_experimental = nan(N,1);
HR_acuerdo_ventanas_original = nan(N,1);
HR_acuerdo_ventanas_experimental = nan(N,1);
HR_dispersion_ventanas_original = nan(N,1);
HR_dispersion_ventanas_experimental = nan(N,1);
HR_RMS_fase_original = nan(N,1);
HR_RMS_fase_experimental = nan(N,1);
HR_fusion_count_original = nan(N,1);
HR_fusion_count_experimental = nan(N,1);
HR_fallback_original = nan(N,1);
HR_fallback_experimental = nan(N,1);
HR_ventanas_bpm_original = repmat({''},N,1);
HR_ventanas_bpm_experimental = repmat({''},N,1);

Estado_original = repmat({''},N,1);
Estado_experimental = repmat({''},N,1);
Error = repmat({''},N,1);

% Versiones independientes
assert(exist('procesar_signos_vitales_original','file') == 2, ...
    'No se encuentra el algoritmo original.');

assert(exist('procesar_signos_vitales_fixed','file') == 2, ...
    'No se encuentra el algoritmo experimental.');

%% Procesamiento
for k = 1:N

    archivo = archivos{k};
    fprintf('\n[%d/%d] Procesando %s\n', k, N, archivo);

    try
        % 1. Cargar captura
        S = load(archivo);

        if isfield(S,'rx_continuous')
            rx = S.rx_continuous;
        elseif isfield(S,'rx')
            rx = S.rx;
        else
            error('No existe rx_continuous ni rx.');
        end

        if ~isfield(S,'cfg')
            error('No existe cfg.');
        end

        cfg = S.cfg;
        cfg.debug_plots = false;

        clear S

        % 2. Frontend comun
        [iq, Fs] = frontend_dsp(rx, cfg);

        clear rx

        % 3. Filtros identicos a replotear.m
        dHP = designfilt('highpassiir', ...
            'FilterOrder',4, ...
            'HalfPowerFrequency',0.04, ...
            'SampleRate',Fs);

        dLP_rr = designfilt('lowpassiir', ...
            'FilterOrder',4, ...
            'HalfPowerFrequency',0.8, ...
            'SampleRate',Fs);

        % 4. Algoritmo original
        try
            [rr1, hr1, info1] = ...
                procesar_signos_vitales_original( ...
                iq, Fs, dHP, dLP_rr);

            RR_original(k) = rr1;
            HR_original(k) = hr1;

            Acuerdo_original(k) = ...
                campo_num(info1,'hr_agreement');

            Candidato_original(k) = ...
                campo_num(info1,'hr_candidate_bpm');

            Confianza_original(k) = ...
                campo_num(info1,'hr_low_confidence');

            % Diagnosticos HR: guardar los estimadores internos aunque bpm_est
            % sea NaN por rechazo de confianza.
            HR_prominencia_original(k) = campo_num(info1,'hr_prominence');
            HR_score_ratio_original(k) = campo_num(info1,'hr_score_ratio');
            HR_espectral_bruto_original(k) = campo_num(info1,'hr_spectral_bpm_raw');
            HR_espectral_fundamental_original(k) = campo_num(info1,'hr_spectral_bpm_fundamental');
            HR_ratio_subarmonico_original(k) = campo_num(info1,'hr_subharmonic_ratio');
            HR_autocorr_bpm_original(k) = campo_num(info1,'hr_autocorr_bpm');
            HR_autocorr_pico_original(k) = campo_num(info1,'hr_autocorr_peak');
            HR_consenso_ventanas_original(k) = campo_num(info1,'hr_folded_bpm');
            HR_acuerdo_ventanas_original(k) = campo_num(info1,'hr_folded_agreement');
            HR_dispersion_ventanas_original(k) = campo_num(info1,'hr_folded_spread');
            HR_RMS_fase_original(k) = campo_num(info1,'hr_signal_rms_rad');
            HR_fusion_count_original(k) = campo_num(info1,'hr_fusion_count');
            HR_fallback_original(k) = campo_num(info1,'hr_fallback');
            HR_ventanas_bpm_original{k} = campo_lista(info1,'hr_window_bpm');

            Estado_original{k} = estado_hr(hr1);

        catch ME1
            Estado_original{k} = 'ERROR';
            Error{k} = ['Original: ' ME1.message];
        end

        % 5. Algoritmo experimental
        try
            [rr2, hr2, info2] = ...
                procesar_signos_vitales_fixed( ...
                iq, Fs, dHP, dLP_rr);

            RR_experimental(k) = rr2;
            HR_experimental(k) = hr2;

            Acuerdo_experimental(k) = ...
                campo_num(info2,'hr_agreement');

            Candidato_experimental(k) = ...
                campo_num(info2,'hr_candidate_bpm');

            Confianza_experimental(k) = ...
                campo_num(info2,'hr_low_confidence');

            % Diagnosticos HR: guardar los estimadores internos aunque bpm_est
            % sea NaN por rechazo de confianza.
            HR_prominencia_experimental(k) = campo_num(info2,'hr_prominence');
            HR_score_ratio_experimental(k) = campo_num(info2,'hr_score_ratio');
            HR_espectral_bruto_experimental(k) = campo_num(info2,'hr_spectral_bpm_raw');
            HR_espectral_fundamental_experimental(k) = campo_num(info2,'hr_spectral_bpm_fundamental');
            HR_ratio_subarmonico_experimental(k) = campo_num(info2,'hr_subharmonic_ratio');
            HR_autocorr_bpm_experimental(k) = campo_num(info2,'hr_autocorr_bpm');
            HR_autocorr_pico_experimental(k) = campo_num(info2,'hr_autocorr_peak');
            HR_consenso_ventanas_experimental(k) = campo_num(info2,'hr_folded_bpm');
            HR_acuerdo_ventanas_experimental(k) = campo_num(info2,'hr_folded_agreement');
            HR_dispersion_ventanas_experimental(k) = campo_num(info2,'hr_folded_spread');
            HR_RMS_fase_experimental(k) = campo_num(info2,'hr_signal_rms_rad');
            HR_fusion_count_experimental(k) = campo_num(info2,'hr_fusion_count');
            HR_fallback_experimental(k) = campo_num(info2,'hr_fallback');
            HR_ventanas_bpm_experimental{k} = campo_lista(info2,'hr_window_bpm');

            Estado_experimental{k} = estado_hr(hr2);

        catch ME2
            Estado_experimental{k} = 'ERROR';

            if isempty(Error{k})
                Error{k} = ['Experimental: ' ME2.message];
            else
                Error{k} = [Error{k} ...
                    ' | Experimental: ' ME2.message];
            end
        end

        clear iq info1 info2 dHP dLP_rr

    catch ME
        Estado_original{k} = 'ERROR';
        Estado_experimental{k} = 'ERROR';
        Error{k} = ME.message;
    end

    fprintf('Original: RR %.2f | HR %.2f\n', ...
        RR_original(k), HR_original(k));

    fprintf('Experimental: RR %.2f | HR %.2f\n', ...
        RR_experimental(k), HR_experimental(k));

    % Guardado de respaldo en cada iteracion
    save('comparacion_ECHO_diagnosticos_progreso.mat', ...
        'archivos','sujetos','muestras', ...
        'RR_original','HR_original', ...
        'RR_experimental','HR_experimental', ...
        'Acuerdo_original','Acuerdo_experimental', ...
        'Candidato_original','Candidato_experimental', ...
        'Confianza_original','Confianza_experimental', ...
        'HR_prominencia_original','HR_prominencia_experimental', ...
        'HR_score_ratio_original','HR_score_ratio_experimental', ...
        'HR_espectral_bruto_original','HR_espectral_bruto_experimental', ...
        'HR_espectral_fundamental_original','HR_espectral_fundamental_experimental', ...
        'HR_ratio_subarmonico_original','HR_ratio_subarmonico_experimental', ...
        'HR_autocorr_bpm_original','HR_autocorr_bpm_experimental', ...
        'HR_autocorr_pico_original','HR_autocorr_pico_experimental', ...
        'HR_consenso_ventanas_original','HR_consenso_ventanas_experimental', ...
        'HR_acuerdo_ventanas_original','HR_acuerdo_ventanas_experimental', ...
        'HR_dispersion_ventanas_original','HR_dispersion_ventanas_experimental', ...
        'HR_RMS_fase_original','HR_RMS_fase_experimental', ...
        'HR_fusion_count_original','HR_fusion_count_experimental', ...
        'HR_fallback_original','HR_fallback_experimental', ...
        'HR_ventanas_bpm_original','HR_ventanas_bpm_experimental', ...
        'Estado_original','Estado_experimental', ...
        'Error','k');
end

%% Exportar resultados
Resultados = table( ...
    archivos, sujetos, muestras, ...
    RR_original, HR_original, ...
    RR_experimental, HR_experimental, ...
    Acuerdo_original, Acuerdo_experimental, ...
    Candidato_original, Candidato_experimental, ...
    Confianza_original, Confianza_experimental, ...
    HR_prominencia_original, HR_prominencia_experimental, ...
    HR_score_ratio_original, HR_score_ratio_experimental, ...
    HR_espectral_bruto_original, HR_espectral_bruto_experimental, ...
    HR_espectral_fundamental_original, HR_espectral_fundamental_experimental, ...
    HR_ratio_subarmonico_original, HR_ratio_subarmonico_experimental, ...
    HR_autocorr_bpm_original, HR_autocorr_bpm_experimental, ...
    HR_autocorr_pico_original, HR_autocorr_pico_experimental, ...
    HR_consenso_ventanas_original, HR_consenso_ventanas_experimental, ...
    HR_acuerdo_ventanas_original, HR_acuerdo_ventanas_experimental, ...
    HR_dispersion_ventanas_original, HR_dispersion_ventanas_experimental, ...
    HR_RMS_fase_original, HR_RMS_fase_experimental, ...
    HR_fusion_count_original, HR_fusion_count_experimental, ...
    HR_fallback_original, HR_fallback_experimental, ...
    HR_ventanas_bpm_original, HR_ventanas_bpm_experimental, ...
    Estado_original, Estado_experimental, Error, ...
    'VariableNames', { ...
    'Archivo','Participante','Muestra', ...
    'RR_original','HR_original', ...
    'RR_experimental','HR_experimental', ...
    'Acuerdo_original','Acuerdo_experimental', ...
    'Candidato_original','Candidato_experimental', ...
    'Confianza_original','Confianza_experimental', ...
    'HR_prominencia_original','HR_prominencia_experimental', ...
    'HR_score_ratio_original','HR_score_ratio_experimental', ...
    'HR_espectral_bruto_original','HR_espectral_bruto_experimental', ...
    'HR_espectral_fundamental_original','HR_espectral_fundamental_experimental', ...
    'HR_ratio_subarmonico_original','HR_ratio_subarmonico_experimental', ...
    'HR_autocorr_bpm_original','HR_autocorr_bpm_experimental', ...
    'HR_autocorr_pico_original','HR_autocorr_pico_experimental', ...
    'HR_consenso_ventanas_original','HR_consenso_ventanas_experimental', ...
    'HR_acuerdo_ventanas_original','HR_acuerdo_ventanas_experimental', ...
    'HR_dispersion_ventanas_original','HR_dispersion_ventanas_experimental', ...
    'HR_RMS_fase_original','HR_RMS_fase_experimental', ...
    'HR_fusion_count_original','HR_fusion_count_experimental', ...
    'HR_fallback_original','HR_fallback_experimental', ...
    'HR_ventanas_bpm_original','HR_ventanas_bpm_experimental', ...
    'Estado_original','Estado_experimental','Error'});

writetable(Resultados,'comparacion_ECHO_diagnosticos.xlsx');
writetable(Resultados,'comparacion_ECHO_diagnosticos.csv');

fprintf('\n========== RESUMEN ==========\n');
fprintf('Capturas: %d\n', N);

fprintf('HR original aceptadas: %d\n', ...
    sum(isfinite(HR_original)));

fprintf('HR experimental aceptadas: %d\n', ...
    sum(isfinite(HR_experimental)));

fprintf('Errores original: %d\n', ...
    sum(strcmp(Estado_original,'ERROR')));

fprintf('Errores experimental: %d\n', ...
    sum(strcmp(Estado_experimental,'ERROR')));

fprintf('Excel generado: comparacion_ECHO_diagnosticos.xlsx\n');

%% Funciones auxiliares
function valor = campo_num(info, campo)
    valor = NaN;

    if isstruct(info) && isfield(info,campo)
        x = info.(campo);

        if (isnumeric(x) || islogical(x)) && isscalar(x)
            valor = double(x);
        end
    end
end

function estado = estado_hr(hr)
    if isfinite(hr)
        estado = 'ACEPTADA';
    else
        estado = 'RECHAZADA';
    end
end

function txt = campo_lista(info, campo)
    % Convierte vectores de diagnostico en texto exportable a Excel/CSV.
    txt = '';

    if ~isstruct(info) || ~isfield(info,campo)
        return;
    end

    x = info.(campo);
    if ~isnumeric(x) || isempty(x)
        return;
    end

    x = x(:);
    partes = strings(1,numel(x));
    for j = 1:numel(x)
        if isfinite(x(j))
            partes(j) = sprintf('%.2f',x(j));
        else
            partes(j) = 'NaN';
        end
    end
    txt = char(strjoin(partes,', '));
end
