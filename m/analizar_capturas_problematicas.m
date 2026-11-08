
%% ANALIZAR_CAPTURAS_PROBLEMATICAS.M
% Diagnostico integral de capturas problematicas ECHO.
%
% Ejecuta el procesamiento original mediante:
% analizar_capturas_problematicas_base.m
%
% Consolida:
% 1. Referencias y errores HR/RR.
% 2. Diagnostico espectral.
% 3. Soporte subarmonico.
% 4. Consistencia temporal por ventanas.
% 5. Autocorrelacion.
% 6. Fusion de estimadores.
% 7. Criterios de confianza.
% 8. Indicadores I/Q y calidad de captura disponibles.
%
% No modifica el procesador de produccion ni las capturas originales.

clc;

%% ================================================================
%  1. CONFIGURACION
% =================================================================

script_base = 'analizar_capturas_problematicas_base.m';
archivo_base = 'diagnostico_capturas_problematicas.mat';

archivo_csv = 'diagnostico_capturas_problematicas_completo.csv';
archivo_mat = 'diagnostico_capturas_problematicas_completo.mat';
archivo_txt = 'informe_diagnostico_capturas.txt';

% Umbrales tomados de procesar_signos_vitales_fixed.m
umbrales = struct();

umbrales.hr_min_bpm = 55;
umbrales.hr_max_bpm = 160;

umbrales.hr_fusion_tolerancia_bpm = 7;
umbrales.hr_ventana_tolerancia_bpm = 5;

umbrales.hr_acuerdo_espectral_ventanas = 0.75;
umbrales.hr_dispersion_espectral_ventanas = 5;

umbrales.hr_acuerdo_autocorr_ventanas = 0.50;
umbrales.hr_dispersion_autocorr_ventanas = 3;
umbrales.hr_autocorr_fusion = 0.22;

umbrales.hr_soporte_subarmonico = 0.80;
umbrales.hr_rms_min = 0.010;
umbrales.hr_autocorr_min = 0.20;
umbrales.hr_prominencia_min = 1.40;
umbrales.hr_ventanas_minimas = 3;

fprintf('\n============================================\n');
fprintf('       DIAGNOSTICO INTEGRAL ECHO\n');
fprintf('============================================\n');

%% ================================================================
%  2. EJECUTAR PROCESAMIENTO ORIGINAL
% =================================================================

if exist(script_base,'file') == 2

    fprintf('\nEjecutando procesamiento original...\n');
    run(script_base);

elseif exist(archivo_base,'file') == 2

    warning(['No se encontro el script base. Se utilizaran ' ...
        'los resultados previamente guardados.']);

else

    error(['No se encontro el script base ni el archivo MAT. ' ...
        'Verifica los nombres y la carpeta de trabajo.']);

end

if exist(archivo_base,'file') ~= 2
    error('No se genero %s.',archivo_base);
end

%% ================================================================
%  3. CARGAR RESULTADOS
% =================================================================

D = load(archivo_base);

campos_necesarios = {'Info_completa','archivos'};

for k = 1:numel(campos_necesarios)
    if ~isfield(D,campos_necesarios{k})
        error('Falta la variable %s en el archivo MAT.', ...
            campos_necesarios{k});
    end
end

Info_completa = D.Info_completa;
archivos = D.archivos;

if isstring(archivos)
    archivos = cellstr(archivos);
elseif ischar(archivos)
    archivos = cellstr(archivos);
end

archivos = archivos(:);
N = numel(archivos);

if ~iscell(Info_completa) || numel(Info_completa) ~= N
    error('Info_completa no coincide con el numero de capturas.');
end

% Referencias guardadas por el procesamiento original.
% Si no existen, se utilizan las referencias conocidas.

if isfield(D,'HR_ref')
    HR_ref = double(D.HR_ref(:));
else
    HR_ref = [72; 81; 80; 68; 65];
end

if isfield(D,'RR_ref')
    RR_ref = double(D.RR_ref(:));
else
    RR_ref = [19; 11; 12; 9; 12];
end

if numel(HR_ref) ~= N || numel(RR_ref) ~= N
    error('Las referencias no coinciden con el numero de capturas.');
end

if isfield(D,'T')
    T_original = D.T;
else
    T_original = table();
end

%% ================================================================
%  4. INICIALIZAR VARIABLES
% =================================================================

HR_est = nan(N,1);
HR_candidate = nan(N,1);
HR_candidate_error = nan(N,1);
HR_error = nan(N,1);

HR_raw = nan(N,1);
HR_fundamental = nan(N,1);
HR_subharmonic = nan(N,1);
HR_half = nan(N,1);
HR_subharmonic_applied = nan(N,1);

HR_prominence = nan(N,1);
HR_score_ratio = nan(N,1);
HR_autocorr = nan(N,1);
HR_autocorr_peak = nan(N,1);

HR_folded = nan(N,1);
HR_folded_agreement = nan(N,1);
HR_folded_spread = nan(N,1);
HR_signal_rms = nan(N,1);

HR_fusion_count = nan(N,1);
HR_fusion_saved = nan(N,1);
HR_fusion_median = nan(N,1);

HR_path_spectral_ac = nan(N,1);
HR_path_spectral_windows = nan(N,1);
HR_path_ac_windows = nan(N,1);

HR_low_stored = nan(N,1);
HR_low_calculated = nan(N,1);

HR_pass_presence = nan(N,1);
HR_pass_periodicity = nan(N,1);
HR_pass_prominence = nan(N,1);
HR_pass_windows = nan(N,1);
HR_pass_fusion = nan(N,1);

HR_n_windows = nan(N,1);
HR_window_text = cell(N,1);
HR_state = cell(N,1);
HR_reason = cell(N,1);
HR_subharmonic_reason = cell(N,1);

RR_est = nan(N,1);
RR_candidate = nan(N,1);
RR_candidate_error = nan(N,1);
RR_error = nan(N,1);
RR_low = nan(N,1);
RR_prominence = nan(N,1);
RR_score_ratio = nan(N,1);
RR_window_agreement = nan(N,1);
RR_amplitude_ratio = nan(N,1);
RR_consensus_prominence = nan(N,1);
RR_window_text = cell(N,1);
RR_state = cell(N,1);

center_offset = nan(N,1);
arc_span = nan(N,1);
circle_radius = nan(N,1);
carrier_offset = nan(N,1);
lo_drift = nan(N,1);

capture_valid = nan(N,1);
overruns = nan(N,1);
saturation = nan(N,1);
incomplete_blocks = nan(N,1);
total_samples = nan(N,1);
peak_global = nan(N,1);

baseline_peak_off = nan(N,1);
baseline_peak_on = nan(N,1);
baseline_rms_off = nan(N,1);
baseline_rms_on = nan(N,1);

tx_gain_actual = nan(N,1);
rx_gain_actual = nan(N,1);

%% ================================================================
%  5. ANALISIS INDIVIDUAL DE CADA CAPTURA
% =================================================================

for k = 1:N

    V = Info_completa{k};

    if ~isstruct(V)
        warning('La captura %s no contiene una estructura valida.', ...
            archivos{k});
        continue;
    end

    fprintf('\n--------------------------------------------\n');
    fprintf('CAPTURA: %s\n',archivos{k});
    fprintf('--------------------------------------------\n');

    %% 5.1 Referencias y estimaciones

    HR_candidate(k) = getnum(V,'hr_candidate_bpm');
    RR_candidate(k) = getnum(V,'rr_candidate_rpm');

    HR_low_stored(k) = getnum(V,'hr_low_confidence');
    RR_low(k) = getnum(V,'rr_low_confidence');

    if isfinite(HR_low_stored(k)) && ...
            HR_low_stored(k) == 0 && isfinite(HR_candidate(k))
        HR_est(k) = HR_candidate(k);
    end

    if isfinite(RR_low(k)) && ...
            RR_low(k) == 0 && isfinite(RR_candidate(k))
        RR_est(k) = RR_candidate(k);
    end

    if isfinite(HR_est(k)) && isfinite(HR_ref(k))
        HR_error(k) = abs(HR_est(k)-HR_ref(k));
    end

    if isfinite(HR_candidate(k)) && isfinite(HR_ref(k))
        HR_candidate_error(k) = ...
            abs(HR_candidate(k)-HR_ref(k));
    end

    if isfinite(RR_est(k)) && isfinite(RR_ref(k))
        RR_error(k) = abs(RR_est(k)-RR_ref(k));
    end

    if isfinite(RR_candidate(k)) && isfinite(RR_ref(k))
        RR_candidate_error(k) = ...
            abs(RR_candidate(k)-RR_ref(k));
    end

    %% 5.2 Diagnostico espectral

    HR_raw(k) = getnum(V,'hr_spectral_bpm_raw');
    HR_fundamental(k) = ...
        getnum(V,'hr_spectral_bpm_fundamental');

    HR_subharmonic(k) = ...
        getnum(V,'hr_subharmonic_ratio');

    HR_half(k) = HR_raw(k)/2;

    HR_prominence(k) = getnum(V,'hr_prominence');
    HR_score_ratio(k) = getnum(V,'hr_score_ratio');

    if isfinite(HR_raw(k)) && isfinite(HR_subharmonic(k))

        if HR_subharmonic(k) < ...
                umbrales.hr_soporte_subarmonico

            HR_subharmonic_applied(k) = 0;
            HR_subharmonic_reason{k} = ...
                'Soporte subarmonico inferior a 0.80';

        elseif HR_half(k) < umbrales.hr_min_bpm

            HR_subharmonic_applied(k) = 0;
            HR_subharmonic_reason{k} = ...
                'La mitad del candidato queda debajo de 55 bpm';

        else

            HR_subharmonic_applied(k) = 1;
            HR_subharmonic_reason{k} = ...
                'Se cumplen las condiciones de correccion';

        end

    else

        HR_subharmonic_reason{k} = 'Datos no disponibles';

    end

    %% 5.3 Diagnostico temporal

    HR_window = getvec(V,'hr_window_bpm');
    RR_window = getvec(V,'rr_window_rpm');

    HR_n_windows(k) = numel(HR_window);

    HR_window_text{k} = vectext(HR_window);
    RR_window_text{k} = vectext(RR_window);

    HR_folded(k) = getnum(V,'hr_folded_bpm');
    HR_folded_agreement(k) = ...
        getnum(V,'hr_folded_agreement');

    HR_folded_spread(k) = ...
        getnum(V,'hr_folded_spread');

    %% 5.4 Autocorrelacion y amplitud de fase

    HR_autocorr(k) = getnum(V,'hr_autocorr_bpm');
    HR_autocorr_peak(k) = getnum(V,'hr_autocorr_peak');

    HR_signal_rms(k) = getnum(V,'hr_signal_rms_rad');

    %% 5.5 Reconstruccion de las tres vias de fusion

    espectral = HR_fundamental(k);
    autocorr = HR_autocorr(k);
    temporal = HR_folded(k);

    acuerdo = HR_folded_agreement(k);
    dispersion = HR_folded_spread(k);
    pico_ac = HR_autocorr_peak(k);

    % Via 1: espectral + autocorrelacion

    if isfinite(espectral) && isfinite(autocorr)

        HR_path_spectral_ac(k) = double( ...
            abs(espectral-autocorr) <= ...
            umbrales.hr_fusion_tolerancia_bpm);

    end

    % Via 2: espectral + ventanas

    if all(isfinite([espectral temporal acuerdo dispersion]))

        HR_path_spectral_windows(k) = double( ...
            abs(espectral-temporal) <= ...
                umbrales.hr_fusion_tolerancia_bpm && ...
            acuerdo >= umbrales.hr_acuerdo_espectral_ventanas && ...
            dispersion <= umbrales.hr_dispersion_espectral_ventanas);

    end

    % Via 3: autocorrelacion + ventanas

    if all(isfinite([autocorr temporal acuerdo dispersion pico_ac]))

        HR_path_ac_windows(k) = double( ...
            abs(autocorr-temporal) <= ...
                umbrales.hr_fusion_tolerancia_bpm && ...
            acuerdo >= umbrales.hr_acuerdo_autocorr_ventanas && ...
            dispersion <= umbrales.hr_dispersion_autocorr_ventanas && ...
            pico_ac >= umbrales.hr_autocorr_fusion);

    end

    % Reconstruccion exacta de la union de estimadores.

    usar = false(1,3);

    if HR_path_spectral_ac(k) == 1
        usar([1 2]) = true;
    end

    if HR_path_spectral_windows(k) == 1
        usar([1 3]) = true;
    end

    if HR_path_ac_windows(k) == 1
        usar([2 3]) = true;
    end

  estimadores = [espectral autocorr temporal];

usar = usar & isfinite(estimadores);

if any(isfinite([HR_path_spectral_ac(k), ...
                 HR_path_spectral_windows(k), ...
                 HR_path_ac_windows(k)]))

    HR_fusion_count(k) = nnz(usar);

    if any(usar)
        HR_fusion_median(k) = median(estimadores(usar));
    end

end

    HR_fusion_saved(k) = getnum(V,'hr_fusion_count');

    %% 5.6 Criterios individuales de confianza

    if isfinite(HR_signal_rms(k))
        HR_pass_presence(k) = double( ...
            HR_signal_rms(k) >= umbrales.hr_rms_min);
    end

    if isfinite(HR_autocorr_peak(k))
        HR_pass_periodicity(k) = double( ...
            HR_autocorr_peak(k) >= umbrales.hr_autocorr_min);
    end

    if isfinite(HR_prominence(k))
        HR_pass_prominence(k) = double( ...
            HR_prominence(k) >= umbrales.hr_prominencia_min);
    end

    if isfinite(HR_n_windows(k))
        HR_pass_windows(k) = double( ...
            HR_n_windows(k) >= umbrales.hr_ventanas_minimas);
    end

    if isfinite(HR_fusion_count(k))
        HR_pass_fusion(k) = double(HR_fusion_count(k) >= 2);
    end

    if all(isfinite([HR_n_windows(k), ...
                     HR_signal_rms(k), ...
                     HR_autocorr_peak(k), ...
                     HR_prominence(k), ...
                     HR_fusion_count(k)]))

        HR_low_calculated(k) = double( ...
            HR_n_windows(k) < umbrales.hr_ventanas_minimas || ...
            HR_signal_rms(k) < umbrales.hr_rms_min || ...
            HR_autocorr_peak(k) < umbrales.hr_autocorr_min || ...
            HR_prominence(k) < umbrales.hr_prominencia_min || ...
            HR_fusion_count(k) < 2);

    end

    %% 5.7 Explicacion automatica del resultado

    motivos = {};

    if HR_pass_windows(k) == 0
        motivos{end+1} = 'Menos de tres ventanas validas';
    end

    if HR_pass_presence(k) == 0
        motivos{end+1} = 'RMS de fase insuficiente';
    end

    if HR_pass_periodicity(k) == 0
        motivos{end+1} = 'Periodicidad insuficiente';
    end

    if HR_pass_prominence(k) == 0
        motivos{end+1} = 'Prominencia espectral insuficiente';
    end

    if HR_pass_fusion(k) == 0
        motivos{end+1} = 'No se obtuvo una fusion valida';
    end

    if isempty(motivos)
        HR_reason{k} = ...
            'Todos los criterios generales de confianza se cumplen';
    else
        HR_reason{k} = strjoin(motivos,'; ');
    end

    if HR_low_stored(k) == 1
        HR_state{k} = 'RECHAZADA';
    elseif isfinite(HR_est(k))
        HR_state{k} = 'ACEPTADA';
    else
        HR_state{k} = 'SIN ESTIMACION';
    end

    %% 5.8 Indicadores RR adicionales

    RR_prominence(k) = getnum(V,'rr_prominence');
    RR_score_ratio(k) = getnum(V,'rr_score_ratio');

    RR_window_agreement(k) = ...
        getnum(V,'rr_window_agreement');

    RR_amplitude_ratio(k) = ...
        getnum(V,'rr_amplitude_ratio');

    RR_consensus_prominence(k) = ...
        getnum(V,'rr_consensus_prominence');

    if RR_low(k) == 1
        RR_state{k} = 'RECHAZADA';
    elseif isfinite(RR_est(k))
        RR_state{k} = 'ACEPTADA';
    else
        RR_state{k} = 'SIN ESTIMACION';
    end

    %% 5.9 Indicadores I/Q

    center_offset(k) = getnum(V,'center_offset_ratio');
    arc_span(k) = getnum(V,'arc_span_rad');
    circle_radius(k) = getnum(V,'circle_radius');
    carrier_offset(k) = getnum(V,'carrier_offset_hz');
    lo_drift(k) = getnum(V,'lo_drift_pp_rad');

    %% 5.10 Integridad de captura, si esta guardada

    capture_valid(k) = getnum(V,'capture_valid');
    overruns(k) = getnum(V,'overruns');
    saturation(k) = getnum(V,'saturation');

    incomplete_blocks(k) = ...
        getnum(V,'n_incomplete_blocks');

    total_samples(k) = getnum(V,'total_actual_samples');
    peak_global(k) = getnum(V,'peak_global');

    baseline_peak_off(k) = getnum(V,'baseline_peak_off');
    baseline_peak_on(k) = getnum(V,'baseline_peak_on');
    baseline_rms_off(k) = getnum(V,'baseline_rms_off');
    baseline_rms_on(k) = getnum(V,'baseline_rms_on');

    tx_gain_actual(k) = getnum(V,'tx_gain_actual');
    rx_gain_actual(k) = getnum(V,'rx_gain_actual');

    %% 5.11 Salida individual

    fprintf('HR referencia: %.2f bpm\n',HR_ref(k));
    fprintf('HR candidato: %.2f bpm\n',HR_candidate(k));
    fprintf('HR final: %.2f bpm\n',HR_est(k));
    fprintf('Error HR aceptado: %.2f bpm\n',HR_error(k));

    fprintf('HR espectral bruto/fundamental: %.2f / %.2f\n', ...
        HR_raw(k),HR_fundamental(k));

    fprintf('Soporte subarmonico: %.4f\n',HR_subharmonic(k));
    fprintf('Prominencia: %.4f\n',HR_prominence(k));
    fprintf('Score ratio: %.4f\n',HR_score_ratio(k));

    fprintf('Autocorrelacion: %.2f bpm (pico %.4f)\n', ...
        HR_autocorr(k),HR_autocorr_peak(k));

    fprintf('Consenso temporal: %.2f bpm\n',HR_folded(k));
    fprintf('Acuerdo: %.2f\n',HR_folded_agreement(k));
    fprintf('Dispersion: %.2f bpm\n',HR_folded_spread(k));

    fprintf('Ventanas: %s\n',HR_window_text{k});

    fprintf('Fusion calculada/guardada: %.0f / %.0f\n', ...
        HR_fusion_count(k),HR_fusion_saved(k));

    fprintf('Estado HR: %s\n',HR_state{k});
    fprintf('Diagnostico: %s\n',HR_reason{k});

    fprintf('RR referencia/estimada: %.2f / %.2f rpm\n', ...
        RR_ref(k),RR_est(k));

end

%% ================================================================
%  6. CONSTRUIR TABLA CONSOLIDADA
% =================================================================

Tdiag = table( ...
    archivos, ...
    HR_ref, HR_est, HR_error, ...
    HR_candidate, HR_candidate_error, HR_state, ...
    HR_low_stored, HR_low_calculated, ...
    HR_raw, HR_fundamental, ...
    HR_subharmonic, HR_half, ...
    HR_subharmonic_applied, HR_subharmonic_reason, ...
    HR_prominence, HR_score_ratio, ...
    HR_autocorr, HR_autocorr_peak, ...
    HR_folded, HR_folded_agreement, HR_folded_spread, ...
    HR_window_text, HR_n_windows, HR_signal_rms, ...
    HR_path_spectral_ac, ...
    HR_path_spectral_windows, ...
    HR_path_ac_windows, ...
    HR_fusion_count, HR_fusion_saved, HR_fusion_median, ...
    HR_pass_presence, HR_pass_periodicity, ...
    HR_pass_prominence, HR_pass_windows, HR_pass_fusion, ...
    HR_reason, ...
    RR_ref, RR_est, RR_error, ...
    RR_candidate, RR_candidate_error, RR_state, RR_low, ...
    RR_prominence, RR_score_ratio, ...
    RR_window_text, RR_window_agreement, ...
    RR_amplitude_ratio, RR_consensus_prominence, ...
    center_offset, arc_span, circle_radius, ...
    carrier_offset, lo_drift, ...
    capture_valid, overruns, saturation, ...
    incomplete_blocks, total_samples, peak_global, ...
    baseline_peak_off, baseline_peak_on, ...
    baseline_rms_off, baseline_rms_on, ...
    tx_gain_actual, rx_gain_actual, ...
    'VariableNames',{ ...
    'Archivo', ...
    'HR_ref_bpm','HR_est_bpm','HR_error_abs_bpm', ...
    'HR_candidate_bpm','HR_candidate_error_abs_bpm','HR_estado', ...
    'HR_low_confidence','HR_low_confidence_calculada', ...
    'HR_spectral_raw_bpm','HR_spectral_fundamental_bpm', ...
    'HR_subharmonic_ratio','HR_half_bpm', ...
    'HR_subharmonic_applied','HR_subharmonic_reason', ...
    'HR_prominence','HR_score_ratio', ...
    'HR_autocorr_bpm','HR_autocorr_peak', ...
    'HR_folded_bpm','HR_folded_agreement','HR_folded_spread', ...
    'HR_window_bpm','HR_n_windows','HR_signal_rms_rad', ...
    'Fusion_spectral_autocorr', ...
    'Fusion_spectral_windows', ...
    'Fusion_autocorr_windows', ...
    'HR_fusion_count_calculated','HR_fusion_count_saved', ...
    'HR_fusion_median_calculated', ...
    'HR_pass_presence','HR_pass_periodicity', ...
    'HR_pass_prominence','HR_pass_windows','HR_pass_fusion', ...
    'HR_diagnosis', ...
    'RR_ref_rpm','RR_est_rpm','RR_error_abs_rpm', ...
    'RR_candidate_rpm','RR_candidate_error_abs_rpm', ...
    'RR_estado','RR_low_confidence', ...
    'RR_prominence','RR_score_ratio', ...
    'RR_window_rpm','RR_window_agreement', ...
    'RR_amplitude_ratio','RR_consensus_prominence', ...
    'center_offset_ratio','arc_span_rad','circle_radius', ...
    'carrier_offset_hz','lo_drift_pp_rad', ...
    'capture_valid','overruns','saturation', ...
    'n_incomplete_blocks','total_actual_samples','peak_global', ...
    'baseline_peak_off','baseline_peak_on', ...
    'baseline_rms_off','baseline_rms_on', ...
    'tx_gain_actual','rx_gain_actual'});

%% ================================================================
%  7. EXPORTAR CSV
% =================================================================

writetable(Tdiag,archivo_csv);

fprintf('\nCSV generado: %s\n',archivo_csv);

%% ================================================================
%  8. EXPORTAR INFORME TXT
% =================================================================

fid = fopen(archivo_txt,'w');

if fid < 0
    error('No se pudo crear el informe TXT.');
end

fprintf(fid,'DIAGNOSTICO INTEGRAL DE CAPTURAS ECHO\n');
fprintf(fid,'=====================================\n\n');

fprintf(fid,'Numero de capturas: %d\n',N);
fprintf(fid,'Fecha: %s\n\n',datestr(now));

fprintf(fid,'UMBRALES UTILIZADOS\n');
fprintf(fid,'HR banda: %.1f - %.1f bpm\n', ...
    umbrales.hr_min_bpm,umbrales.hr_max_bpm);
fprintf(fid,'Tolerancia fusion: %.1f bpm\n', ...
    umbrales.hr_fusion_tolerancia_bpm);
fprintf(fid,'Soporte subarmonico: %.2f\n', ...
    umbrales.hr_soporte_subarmonico);
fprintf(fid,'RMS minimo: %.3f rad\n',umbrales.hr_rms_min);
fprintf(fid,'Autocorrelacion minima: %.2f\n', ...
    umbrales.hr_autocorr_min);
fprintf(fid,'Prominencia minima: %.2f\n\n', ...
    umbrales.hr_prominencia_min);

for k = 1:N

    fprintf(fid,'\n=====================================\n');
    fprintf(fid,'CAPTURA: %s\n',archivos{k});
    fprintf(fid,'=====================================\n');

    fprintf(fid,'HR referencia: %.3f bpm\n',HR_ref(k));
    fprintf(fid,'HR estimada: %.3f bpm\n',HR_est(k));
    fprintf(fid,'Error HR aceptado: %.3f bpm\n',HR_error(k));
    fprintf(fid,'HR candidato: %.3f bpm\n',HR_candidate(k));
    fprintf(fid,'Error candidato: %.3f bpm\n',HR_candidate_error(k));
    fprintf(fid,'Estado: %s\n\n',HR_state{k});

    fprintf(fid,'ESPECTRO HR\n');
    fprintf(fid,'Bruto: %.4f bpm\n',HR_raw(k));
    fprintf(fid,'Fundamental: %.4f bpm\n',HR_fundamental(k));
    fprintf(fid,'Soporte subarmonico: %.4f\n',HR_subharmonic(k));
    fprintf(fid,'Mitad del candidato: %.4f bpm\n',HR_half(k));
    fprintf(fid,'Correccion aplicada: %.0f\n', ...
        HR_subharmonic_applied(k));
    fprintf(fid,'Motivo: %s\n',HR_subharmonic_reason{k});
    fprintf(fid,'Prominencia: %.4f\n',HR_prominence(k));
    fprintf(fid,'Score ratio: %.4f\n\n',HR_score_ratio(k));

    fprintf(fid,'CONSISTENCIA TEMPORAL\n');
    fprintf(fid,'Ventanas: %s\n',HR_window_text{k});
    fprintf(fid,'Consenso: %.4f bpm\n',HR_folded(k));
    fprintf(fid,'Acuerdo: %.4f\n',HR_folded_agreement(k));
    fprintf(fid,'Dispersion: %.4f bpm\n\n',HR_folded_spread(k));

    fprintf(fid,'AUTOCORRELACION Y FASE\n');
    fprintf(fid,'Autocorrelacion: %.4f bpm\n',HR_autocorr(k));
    fprintf(fid,'Pico autocorrelacion: %.4f\n',HR_autocorr_peak(k));
    fprintf(fid,'RMS fase: %.5f rad\n\n',HR_signal_rms(k));

    fprintf(fid,'VIAS DE FUSION\n');
    fprintf(fid,'Espectral + autocorrelacion: %.0f\n', ...
        HR_path_spectral_ac(k));
    fprintf(fid,'Espectral + ventanas: %.0f\n', ...
        HR_path_spectral_windows(k));
    fprintf(fid,'Autocorrelacion + ventanas: %.0f\n', ...
        HR_path_ac_windows(k));
    fprintf(fid,'Fusion calculada: %.0f\n',HR_fusion_count(k));
    fprintf(fid,'Fusion guardada: %.0f\n',HR_fusion_saved(k));
    fprintf(fid,'Mediana fusion calculada: %.4f bpm\n\n', ...
        HR_fusion_median(k));

    fprintf(fid,'CRITERIOS DE CONFIANZA\n');
    fprintf(fid,'Presencia: %.0f\n',HR_pass_presence(k));
    fprintf(fid,'Periodicidad: %.0f\n',HR_pass_periodicity(k));
    fprintf(fid,'Prominencia: %.0f\n',HR_pass_prominence(k));
    fprintf(fid,'Ventanas: %.0f\n',HR_pass_windows(k));
    fprintf(fid,'Fusion: %.0f\n',HR_pass_fusion(k));
    fprintf(fid,'Baja confianza guardada: %.0f\n',HR_low_stored(k));
    fprintf(fid,'Baja confianza reconstruida: %.0f\n', ...
        HR_low_calculated(k));
    fprintf(fid,'Diagnostico: %s\n\n',HR_reason{k});

    fprintf(fid,'RESPIRACION\n');
    fprintf(fid,'RR referencia: %.3f rpm\n',RR_ref(k));
    fprintf(fid,'RR estimada: %.3f rpm\n',RR_est(k));
    fprintf(fid,'Error RR: %.3f rpm\n',RR_error(k));
    fprintf(fid,'Ventanas RR: %s\n',RR_window_text{k});
    fprintf(fid,'Prominencia RR: %.4f\n',RR_prominence(k));
    fprintf(fid,'Score ratio RR: %.4f\n',RR_score_ratio(k));

    fprintf(fid,'\nINTEGRIDAD DE CAPTURA\n');
    fprintf(fid,'Capture valid: %.0f\n',capture_valid(k));
    fprintf(fid,'Overruns: %.0f\n',overruns(k));
    fprintf(fid,'Saturation: %.0f\n',saturation(k));
    fprintf(fid,'Bloques incompletos: %.0f\n',incomplete_blocks(k));
    fprintf(fid,'Muestras: %.0f\n',total_samples(k));
    fprintf(fid,'Peak global: %.5f\n',peak_global(k));

end

fprintf(fid,'\n=====================================\n');
fprintf(fid,'FIN DEL INFORME\n');
fprintf(fid,'=====================================\n');

fclose(fid);

fprintf('Informe TXT generado: %s\n',archivo_txt);

%% ================================================================
%  9. GUARDAR MAT COMPLETO
% =================================================================

save(archivo_mat, ...
    'Tdiag','Info_completa','archivos', ...
    'HR_ref','RR_ref','umbrales','T_original','-v7.3');

fprintf('MAT completo generado: %s\n',archivo_mat);

%% ================================================================
%  10. RESUMEN FINAL
% =================================================================

fprintf('\n============================================\n');
fprintf('          RESUMEN DEL DIAGNOSTICO\n');
fprintf('============================================\n');

fprintf('Capturas analizadas: %d\n',N);
fprintf('HR aceptadas: %d\n',nnz(strcmp(HR_state,'ACEPTADA')));
fprintf('HR rechazadas: %d\n',nnz(strcmp(HR_state,'RECHAZADA')));

fprintf('\nErrores HR de estimaciones aceptadas:\n');
disp(HR_error);

fprintf('Errores RR de estimaciones aceptadas:\n');
disp(RR_error);

fprintf('\nTabla consolidada:\n');
disp(Tdiag(:,{'Archivo','HR_ref_bpm','HR_est_bpm', ...
    'HR_error_abs_bpm','HR_estado', ...
    'HR_fusion_count_calculated', ...
    'HR_subharmonic_ratio','RR_error_abs_rpm'}));

fprintf('\nARCHIVOS GENERADOS:\n');
fprintf('1. %s\n',archivo_csv);
fprintf('2. %s\n',archivo_mat);
fprintf('3. %s\n',archivo_txt);

fprintf('\nDiagnostico terminado.\n');

%% ================================================================
%  FUNCIONES AUXILIARES
% =================================================================

function x = getnum(S,campo)

    x = NaN;

    if isstruct(S) && isfield(S,campo)

        v = S.(campo);

        if (isnumeric(v) || islogical(v)) && isscalar(v)
            x = double(v);
        end

    end

end

function v = getvec(S,campo)

    v = [];

    if isstruct(S) && isfield(S,campo)

        x = S.(campo);

        if isnumeric(x) || islogical(x)
            v = double(x(:)).';
            v = v(isfinite(v));
        end

    end

end

function t = vectext(v)

    if isempty(v)
        t = '[]';
    else
        t = mat2str(v,5);
    end

end
