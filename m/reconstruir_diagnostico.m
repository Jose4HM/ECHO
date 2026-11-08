
%% RECONSTRUCCION DE DIAGNOSTICOS
% Utiliza exclusivamente datos guardados en vital_info.
% Reproduce las condiciones del procesador definitivo.
% No modifica los archivos MAT.

clc;
clear;

archivos = {
    'S2_4.mat'
    'S3_2.mat'
    'S7_9.mat'
    'S7_10.mat'
    'S8_1.mat'
    'S8_4.mat'
};

hr_min = 55;
hr_max = 160;

N = numel(archivos);

Resumen = cell(N,1);

%% RESULTADOS DE PICOS ESPECTRALES

p_archivo = {};
p_orden = [];
p_frecuencia = [];
p_magnitud = [];
p_prominencia = [];

%% PROCESAMIENTO

for n = 1:N

    fprintf('\n====================================\n');
    fprintf('DIAGNOSTICO %s\n',archivos{n});
    fprintf('====================================\n');

    S = load(archivos{n},'vital_info');
    V = S.vital_info;

    Fs = V.Fs;
    x = V.x_hr(:);
    f = V.f_hr_axis(:);
    P = V.P_hr(:);

    f_est = V.f_hr_est;
    ventanas = V.hr_window_bpm(:);

    %% 1. HR ESPECTRAL BRUTO

    HR_raw = 60*f_est;

    %% 2. SOPORTE SUBARMONICO

    f_sub = f_est/2;

    m_sub = abs(f-f_sub) <= 0.035;
    m_full = abs(f-f_est) <= 0.035;

    if any(m_sub) && any(m_full)
        soporte = max(P(m_sub))/(max(P(m_full))+eps);
    else
        soporte = 0;
    end

    HR_fundamental = HR_raw;

    if f_est/2 >= hr_min/60 && soporte >= 0.80
        HR_fundamental = HR_raw/2;
    end

    %% 3. AUTOCORRELACION
    % Reproduccion de candidato_hr_autocorr

    xc = x-mean(x);

    if numel(xc) < round(3*Fs) || sqrt(mean(xc.^2)) <= eps

        HR_ac = NaN;
        pico_ac = 0;

    else

        r = xcorr(xc,'coeff');
        r = r(numel(xc):end);

        i_lo = max(2,round(60/(hr_max)*Fs)+1);
        i_hi = min(numel(r)-1, ...
            round(60/(hr_min)*Fs)+1);

        if i_hi <= i_lo

            HR_ac = NaN;
            pico_ac = 0;

        else

            idx = i_lo:i_hi;

            loc = idx( ...
                r(idx) >= r(idx-1) & ...
                r(idx) >= r(idx+1));

            if isempty(loc)
                [pico_ac,k] = max(r(idx));
                kk = idx(k);
            else
                [pico_ac,k] = max(r(loc));
                kk = loc(k);
            end

            % Refinamiento parabolico

            d = 0;

            if kk > 1 && kk < numel(r)

                den = r(kk-1)-2*r(kk)+r(kk+1);

                if abs(den) > 1e-12
                    d = 0.5*(r(kk-1)-r(kk+1))/den;
                    d = max(-0.5,min(0.5,d));
                end

            end

            lag = (kk-1)+d;
            HR_ac = 60*Fs/lag;

        end
    end

    %% 4. CONSENSO TEMPORAL

    bpm_win = ventanas(isfinite(ventanas));

    opciones = [];

    for k = 1:numel(bpm_win)

        if bpm_win(k) >= hr_min && bpm_win(k) <= hr_max
            opciones(end+1,:) = [bpm_win(k),k];
        end

        mitad = bpm_win(k)/2;

        if mitad >= hr_min && mitad <= hr_max
            opciones(end+1,:) = [mitad,k];
        end

    end

    if isempty(opciones)

        HR_folded = NaN;
        acuerdo = 0;
        dispersion = Inf;

    else

        rejilla = hr_min:0.25:hr_max;
        cuentas = zeros(size(rejilla));

        for j = 1:numel(rejilla)

            votos = false(1,numel(bpm_win));

            for q = 1:size(opciones,1)

                if abs(opciones(q,1)-rejilla(j)) <= 5
                    votos(opciones(q,2)) = true;
                end

            end

            cuentas(j) = nnz(votos);

        end

        [n_votos,j] = max(cuentas);
        centro_previo = rejilla(j);

        elegidos = NaN(size(bpm_win));

        for k = 1:numel(bpm_win)

            vals = bpm_win(k);
            mitad = bpm_win(k)/2;

            if mitad >= hr_min && mitad <= hr_max
                vals = [vals,mitad];
            end

            [dist,q] = min(abs(vals-centro_previo));

            if dist <= 5
                elegidos(k) = vals(q);
            end

        end

        elegidos = elegidos(isfinite(elegidos));

        acuerdo = n_votos/max(numel(bpm_win),1);

        if isempty(elegidos)

            HR_folded = NaN;
            dispersion = Inf;

        else

            HR_folded = median(elegidos);

            dispersion = 1.4826* ...
                median(abs(elegidos-HR_folded));

        end

    end

    %% 5. TRES VIAS DE FUSION

    est = [HR_fundamental,HR_ac,HR_folded];

    prueba1 = isfinite(est(1)) && ...
              isfinite(est(2)) && ...
              abs(est(1)-est(2)) <= 7;

    prueba2 = isfinite(est(1)) && ...
              isfinite(est(3)) && ...
              abs(est(1)-est(3)) <= 7 && ...
              acuerdo >= 0.75 && dispersion <= 5;

    prueba3 = isfinite(est(2)) && ...
              isfinite(est(3)) && ...
              abs(est(2)-est(3)) <= 7 && ...
              acuerdo >= 0.50 && dispersion <= 3 && ...
              pico_ac >= 0.22;

    usar = false(1,3);

    if prueba1
        usar([1 2]) = true;
    end

    if prueba2
        usar([1 3]) = true;
    end

    if prueba3
        usar([2 3]) = true;
    end

    usar = usar & isfinite(est);

    fusion_count = nnz(usar);

    if any(usar)
        HR_fusion = median(est(usar));
    else
        HR_fusion = NaN;
    end

    %% 6. RMS Y CONFIANZA

    x_c = x-median(x);
    RMS_hr = sqrt(mean(x_c.^2));

    presencia = RMS_hr >= 0.010;
    periodicidad = pico_ac >= 0.20;
    evidencia_espectral = V.hr_prominence >= 1.40;

    fusion_valida = isfinite(HR_fusion) && ...
                    fusion_count >= 2;

    baja_confianza = numel(bpm_win) < 3 || ...
                     ~presencia || ...
                     ~periodicidad || ...
                     ~evidencia_espectral || ...
                     ~fusion_valida;

    if fusion_valida
        HR_candidato_reconstruido = HR_fusion;
    else
        HR_candidato_reconstruido = HR_raw;
    end

    %% 7. PRINCIPALES PICOS ESPECTRALES
    % El selector HR utiliza prominencia como puntuacion.
    % f_rr es NaN en el procesamiento definitivo.

    mask = f >= hr_min/60 & f <= hr_max/60;
    idx_band = find(mask);

    idx_peaks = [];

    for ii = 1:numel(idx_band)

        k = idx_band(ii);

        if k > 1 && k < numel(P) && ...
           P(k) >= P(k-1) && P(k) >= P(k+1)

            idx_peaks(end+1) = k;

        end

    end

    if isempty(idx_peaks)
        [~,q] = max(P(idx_band));
        idx_peaks = idx_band(q);
    end

    scores = zeros(size(idx_peaks));

    for ii = 1:numel(idx_peaks)

        k = idx_peaks(ii);
        f0 = f(k);

        nf_mask = mask & ...
            abs(f-f0) <= 0.35 & ...
            abs(f-f0) >= 0.06;

        if nnz(nf_mask) < 8
            nf_mask = mask;
        end

        if any(nf_mask)
            nf = median(P(nf_mask));
        else
            nf = NaN;
        end

        if ~isfinite(nf) || nf <= 0
            nf = median(P(mask))+eps;
        end

        scores(ii) = P(k)/(nf+eps);

    end

    [scores_sorted,orden] = sort(scores,'descend');

    cantidad = min(10,numel(orden));

    fprintf('\n--- DIEZ PRINCIPALES PICOS ---\n');

    fprintf('Orden   HR (bpm)   Magnitud   Prominencia\n');

    for j = 1:cantidad

        k = idx_peaks(orden(j));

        fprintf('%3d    %8.3f   %9.4f   %10.4f\n', ...
            j,60*f(k),P(k),scores_sorted(j));

        p_archivo{end+1,1} = archivos{n};
        p_orden(end+1,1) = j;
        p_frecuencia(end+1,1) = 60*f(k);
        p_magnitud(end+1,1) = P(k);
        p_prominencia(end+1,1) = scores_sorted(j);

    end

    %% 8. RESULTADOS

    fprintf('\n--- DIAGNOSTICO RECONSTRUIDO ---\n');

    fprintf('HR bruto:              %.4f bpm\n',HR_raw);
    fprintf('Soporte subarmonico:   %.4f\n',soporte);
    fprintf('HR fundamental:       %.4f bpm\n',HR_fundamental);
    fprintf('HR autocorrelacion:   %.4f bpm\n',HR_ac);
    fprintf('Pico autocorrelacion: %.4f\n',pico_ac);
    fprintf('HR consenso:          %.4f bpm\n',HR_folded);
    fprintf('Acuerdo:              %.4f\n',acuerdo);
    fprintf('Dispersion:           %.4f bpm\n',dispersion);
    fprintf('RMS fase:             %.5f rad\n',RMS_hr);

    fprintf('\nPrueba 1: %d\n',prueba1);
    fprintf('Prueba 2: %d\n',prueba2);
    fprintf('Prueba 3: %d\n',prueba3);

    fprintf('Fusion count:         %d\n',fusion_count);
    fprintf('HR fusionada:         %.4f bpm\n',HR_fusion);
    fprintf('Baja confianza:       %d\n',baja_confianza);

    fprintf('Candidato guardado:   %.4f bpm\n', ...
        V.hr_candidate_bpm);

    fprintf('Candidato reconstruido: %.4f bpm\n', ...
        HR_candidato_reconstruido);

    fprintf('Baja confianza guardada: %d\n', ...
        V.hr_low_confidence);

    %% Guardar resumen

    R = struct();

    R.Archivo = archivos{n};
    R.HR_raw = HR_raw;
    R.Soporte = soporte;
    R.HR_fundamental = HR_fundamental;
    R.HR_ac = HR_ac;
    R.Pico_AC = pico_ac;
    R.HR_folded = HR_folded;
    R.Acuerdo = acuerdo;
    R.Dispersion = dispersion;
    R.RMS = RMS_hr;
    R.Prueba1 = prueba1;
    R.Prueba2 = prueba2;
    R.Prueba3 = prueba3;
    R.Fusion_count = fusion_count;
    R.HR_fusion = HR_fusion;
    R.Baja_confianza = baja_confianza;
    R.Candidato_guardado = V.hr_candidate_bpm;
    R.Candidato_reconstruido = HR_candidato_reconstruido;

    Resumen{n} = R;

end

%% EXPORTAR RESULTADOS

T = struct2table([Resumen{:}]);

Picos = table(p_archivo,p_orden,p_frecuencia, ...
    p_magnitud,p_prominencia, ...
    'VariableNames',{'Archivo','Orden','Frecuencia_bpm', ...
    'Magnitud','Prominencia'});

writetable(T,'reconstruccion_seis_capturas.csv');
writetable(Picos,'picos_seis_capturas.csv');

fprintf('\n====================================\n');
fprintf('ARCHIVOS GENERADOS CORRECTAMENTE\n');
fprintf('reconstruccion_S2_4_S3_2.csv\n');
fprintf('picos_S2_4_S3_2.csv\n');
fprintf('====================================\n');
