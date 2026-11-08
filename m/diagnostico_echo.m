
%% ECHO - Diagnose multiple MAT files
clear; clc;

% Write the exact filenames here
archivos = {
    'S2_4.mat'
    'S3_2.mat'
    'S7_9.mat'
    'S7_10.mat'
    'S8_10.mat'
};

% Folder containing the MAT files
carpeta = pwd;

for i = 1:numel(archivos)

    archivo = fullfile(carpeta, archivos{i});

    fprintf('\n\n====================================\n');
    fprintf('CAPTURE: %s\n', archivos{i});
    fprintf('====================================\n');

    if ~isfile(archivo)
        fprintf('ERROR: File not found: %s\n', archivo);
        continue;
    end

    S = load(archivo, 'vital_info', 'bpm', 'rpm', 'cfg');

    % Final estimates
    fprintf('\n--- FINAL OUTPUTS ---\n');
    if isfield(S,'bpm'), fprintf('HR: '); disp(S.bpm); end
    if isfield(S,'rpm'), fprintf('RR: '); disp(S.rpm); end

    % Acquisition diagnostics
    if isfield(S,'cfg')
        C = S.cfg;
        fields = {'capture_valid','saturation', ...
            'overruns','n_incomplete_blocks'};

        fprintf('\n--- ACQUISITION ---\n');
        for k = 1:numel(fields)
            f = fields{k};
            if isfield(C,f)
                fprintf('%s: ',f);
                disp(C.(f));
            end
        end
    end

    % Signal-processing diagnostics
    if isfield(S,'vital_info')
        V = S.vital_info;

        fields = {
            'hr_candidate_bpm'
            'hr_window_bpm'
            'hr_prominence'
            'hr_score_ratio'
            'hr_agreement'
            'hr_low_confidence'
            'rr_candidate_rpm'
            'rr_window_rpm'
            'rr_prominence'
            'rr_score_ratio'
            'rr_window_agreement'
            'rr_low_confidence'
            'center_offset_ratio'
            'arc_span_rad'
            'circle_radius'
            'carrier_offset_hz'
            'lo_drift_pp_rad'
        };

        fprintf('\n--- DIAGNOSTICS ---\n');
        for k = 1:numel(fields)
            f = fields{k};
            if isfield(V,f)
                fprintf('%s: ',f);
                disp(V.(f));
            end
        end
    end
end

fprintf('\n\nALL DIAGNOSTICS COMPLETED.\n');
