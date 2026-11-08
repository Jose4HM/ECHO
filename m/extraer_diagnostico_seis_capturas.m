
%% EXTRAER DIAGNOSTICO DE SEIS CAPTURAS
% No modifica los archivos originales.
% Extrae indicadores guardados en vital_info.

clc;

archivos = {
    'S2_4.mat'
    'S3_2.mat'
    'S7_9.mat'
    'S7_10.mat'
    'S8_1.mat'
    'S8_4.mat'
};

N = numel(archivos);

% Inicializacion
campos = {
    'hr_candidate_bpm'
    'hr_spectral_bpm_raw'
    'hr_spectral_bpm_fundamental'
    'hr_subharmonic_ratio'
    'hr_prominence'
    'hr_score_ratio'
    'hr_autocorr_bpm'
    'hr_autocorr_peak'
    'hr_folded_bpm'
    'hr_folded_agreement'
    'hr_folded_spread'
    'hr_signal_rms_rad'
    'hr_fusion_count'
    'hr_low_confidence'
    'rr_candidate_rpm'
    'rr_prominence'
    'rr_score_ratio'
    'rr_low_confidence'
};

for j = 1:numel(campos)
    datos.(campos{j}) = nan(N,1);
end

HR_ventanas = cell(N,1);
RR_ventanas = cell(N,1);
Picos_espectrales = cell(N,1);

Prueba_1 = nan(N,1);
Prueba_2 = nan(N,1);
Prueba_3 = nan(N,1);

Pasa_presencia = nan(N,1);
Pasa_periodicidad = nan(N,1);
Pasa_prominencia = nan(N,1);
Pasa_ventanas = nan(N,1);
Pasa_fusion = nan(N,1);

Correccion_subarmonica = nan(N,1);

%% Procesamiento individual

for k = 1:N

    fprintf('\nAnalizando %s...\n', archivos{k});

    S = load(archivos{k},'vital_info');

    if ~isfield(S,'vital_info')
        warning('No se encontro vital_info en %s',archivos{k});
        continue;
    end

    V = S.vital_info;

    % Extraer indicadores numericos
    for j = 1:numel(campos)
        nombre = campos{j};

        if isfield(V,nombre)
            valor = V.(nombre);

            if isnumeric(valor) || islogical(valor)
                if isscalar(valor)
                    datos.(nombre)(k) = double(valor);
                end
            end
        end
    end

    % Ventanas HR y RR
    HR_ventanas{k} = vector_texto(V,'hr_window_bpm');
    RR_ventanas{k} = vector_texto(V,'rr_window_rpm');

    %% Evaluar las tres vias de fusion

    esp = getnum(V,'hr_spectral_bpm_fundamental');
    ac = getnum(V,'hr_autocorr_bpm');
    temporal = getnum(V,'hr_folded_bpm');

    acuerdo = getnum(V,'hr_folded_agreement');
    dispersion = getnum(V,'hr_folded_spread');
    pico_ac = getnum(V,'hr_autocorr_peak');

    if all(isfinite([esp ac]))
        Prueba_1(k) = abs(esp-ac) <= 7;
    end

    if all(isfinite([esp temporal acuerdo dispersion]))
        Prueba_2(k) = ...
            abs(esp-temporal) <= 7 && ...
            acuerdo >= 0.75 && dispersion <= 5;
    end

    if all(isfinite([ac temporal acuerdo dispersion pico_ac]))
        Prueba_3(k) = ...
            abs(ac-temporal) <= 7 && ...
            acuerdo >= 0.50 && dispersion <= 3 && ...
            pico_ac >= 0.22;
    end

    %% Criterios generales de confianza

    rms_hr = getnum(V,'hr_signal_rms_rad');
    prom = getnum(V,'hr_prominence');
    ventanas = getvec(V,'hr_window_bpm');
    fusion = getnum(V,'hr_fusion_count');

    if isfinite(rms_hr)
        Pasa_presencia(k) = rms_hr >= 0.010;
    end

    if isfinite(pico_ac)
        Pasa_periodicidad(k) = pico_ac >= 0.20;
    end

    if isfinite(prom)
        Pasa_prominencia(k) = prom >= 1.40;
    end

    if isfield(V,'hr_window_bpm')
        Pasa_ventanas(k) = numel(ventanas) >= 3;
    end

    if isfinite(fusion)
        Pasa_fusion(k) = fusion >= 2;
    end

    %% Correccion subarmonica

    bruto = getnum(V,'hr_spectral_bpm_raw');
    soporte = getnum(V,'hr_subharmonic_ratio');

    if isfinite(bruto) && isfinite(soporte)
        Correccion_subarmonica(k) = ...
            (bruto/2 >= 55) && (soporte >= 0.80);
    end

    %% Diez principales picos espectrales

    Picos_espectrales{k} = obtener_picos(V);

end

%% Construir tabla

Archivo = archivos;

T = table(Archivo, ...
    datos.hr_candidate_bpm, ...
    datos.hr_spectral_bpm_raw, ...
    datos.hr_spectral_bpm_fundamental, ...
    datos.hr_subharmonic_ratio, ...
    Correccion_subarmonica, ...
    datos.hr_prominence, ...
    datos.hr_score_ratio, ...
    datos.hr_autocorr_bpm, ...
    datos.hr_autocorr_peak, ...
    datos.hr_folded_bpm, ...
    datos.hr_folded_agreement, ...
    datos.hr_folded_spread, ...
    HR_ventanas, ...
    datos.hr_signal_rms_rad, ...
    Prueba_1, Prueba_2, Prueba_3, ...
    Pasa_presencia, Pasa_periodicidad, ...
    Pasa_prominencia, Pasa_ventanas, Pasa_fusion, ...
    datos.hr_fusion_count, ...
    datos.hr_low_confidence, ...
    Picos_espectrales, ...
    RR_ventanas, ...
    datos.rr_candidate_rpm, ...
    datos.rr_prominence, ...
    datos.rr_score_ratio, ...
    datos.rr_low_confidence);

T.Properties.VariableNames = {
    'Archivo'
    'HR_candidato'
    'HR_espectral_bruto'
    'HR_fundamental'
    'Soporte_subarmonico'
    'Correccion_subarmonica'
    'Prominencia_HR'
    'Score_ratio_HR'
    'HR_autocorrelacion'
    'Pico_autocorrelacion'
    'HR_consenso_temporal'
    'Acuerdo_temporal'
    'Dispersion_temporal'
    'HR_por_ventanas'
    'RMS_fase_HR'
    'Prueba_1_espectral_AC'
    'Prueba_2_espectral_ventanas'
    'Prueba_3_AC_ventanas'
    'Pasa_presencia'
    'Pasa_periodicidad'
    'Pasa_prominencia'
    'Pasa_numero_ventanas'
    'Pasa_fusion'
    'Fusion_count'
    'Baja_confianza'
    'Diez_picos_espectrales'
    'RR_por_ventanas'
    'RR_candidato'
    'RR_prominencia'
    'RR_score_ratio'
    'RR_baja_confianza'
};

%% Exportar

writetable(T,'diagnostico_seis_capturas.csv');

disp(T);

fprintf('\nCSV generado correctamente:\n');
fprintf('diagnostico_seis_capturas.csv\n');

%% ==========================================================
% FUNCIONES AUXILIARES
% ==========================================================

function x = getnum(S,nombre)

x = NaN;

if isfield(S,nombre)
    v = S.(nombre);

    if (isnumeric(v) || islogical(v)) && isscalar(v)
        x = double(v);
    end
end

end

function v = getvec(S,nombre)

v = [];

if isfield(S,nombre)
    x = S.(nombre);

    if isnumeric(x)
        v = double(x(:)).';
    end
end

end

function t = vector_texto(S,nombre)

v = getvec(S,nombre);

if isempty(v)
    t = 'No disponible';
else
    t = mat2str(v,5);
end

end

function t = obtener_picos(S)

t = 'No disponible';

if ~isfield(S,'P_hr') || ~isfield(S,'f_hr_axis')
    return;
end

P = double(S.P_hr(:));
f = double(S.f_hr_axis(:));

if numel(P) ~= numel(f) || numel(P) < 3
    return;
end

mask = f >= 55/60 & f <= 160/60;
idx = find(mask);

if numel(idx) < 3
    return;
end

% Identificar maximos locales
loc = idx(idx > 1 & idx < numel(P));
loc = loc(P(loc) >= P(loc-1) & P(loc) > P(loc+1));

if isempty(loc)
    [~,q] = max(P(idx));
    loc = idx(q);
end

[~,orden] = sort(P(loc),'descend');
loc = loc(orden(1:min(10,numel(orden))));

partes = cell(numel(loc),1);

for j = 1:numel(loc)
    partes{j} = sprintf('%.3f bpm (%.5g)', ...
        60*f(loc(j)),P(loc(j)));
end

t = strjoin(partes,'; ');

end
