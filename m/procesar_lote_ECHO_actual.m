
%% ECHO - PROCESAMIENTO MASIVO DE CAPTURAS IQ
% Utiliza el frontend y procesador actuales.
% No compara algoritmos ni modifica capturas originales.
% Exporta resultados a Excel, CSV y MAT.
% Este script obtendrá los resultados de todos los mats de la carpeta, usa
% el script replotear.m
clc;
fprintf('\n============================================\n');
fprintf('       ECHO - PROCESAMIENTO POR LOTES\n');
fprintf('============================================\n');

%% 1. VERIFICAR DEPENDENCIAS

assert(exist('replotear','file') == 2, ...
    'No se encuentra replotear.m');

assert(exist('frontend_dsp','file') == 2, ...
    'No se encuentra frontend_dsp.m');

assert(exist('procesar_signos_vitales_fixed','file') == 2, ...
    'No se encuentra procesar_signos_vitales_fixed.m');

assert(nargin('replotear') >= 2, ...
    'replotear debe aceptar el argumento mostrar_figura.');

fprintf('\nFrontend: %s\n', which('frontend_dsp'));
fprintf('Procesador: %s\n', which('procesar_signos_vitales_fixed'));

%% 2. SELECCIONAR CARPETA DE CAPTURAS

captureFolder = uigetdir(pwd, ...
    'Seleccionar carpeta con las capturas IQ originales');

if isequal(captureFolder,0)
    fprintf('Operacion cancelada.\n');
    return;
end

% Buscar MAT en la carpeta y subcarpetas.
files = dir(fullfile(captureFolder,'**','*.mat'));

if isempty(files)
    error('No se encontraron archivos MAT.');
end

%% 3. IDENTIFICAR CAPTURAS VALIDAS

validFiles = {};
skippedFiles = {};
skippedReasons = {};

fprintf('\nRevisando archivos MAT...\n');

for k = 1:numel(files)

    fullName = fullfile(files(k).folder,files(k).name);

    try
        vars = whos('-file',fullName);
        names = {vars.name};

        hasIQ = any(strcmp(names,'rx_continuous')) || ...
                any(strcmp(names,'rx'));

        hasCFG = any(strcmp(names,'cfg'));

        if hasIQ && hasCFG
            validFiles{end+1,1} = fullName;
        else
            skippedFiles{end+1,1} = fullName;
            skippedReasons{end+1,1} = ...
                'No contiene IQ y cfg requeridos';
        end

    catch ME
        skippedFiles{end+1,1} = fullName;
        skippedReasons{end+1,1} = ME.message;
    end
end

N = numel(validFiles);

fprintf('Capturas validas encontradas: %d\n',N);
fprintf('Archivos omitidos: %d\n',numel(skippedFiles));

if N == 0
    error('No se encontraron capturas IQ validas.');
end

%% 4. CREAR CARPETA DE RESULTADOS

timestamp = datestr(now,'yyyymmdd_HHMMSS');

outputFolder = fullfile(captureFolder, ...
    ['ECHO_resultados_' timestamp]);

mkdir(outputFolder);

%% 5. PREPARAR RESULTADOS

results = struct( ...
    'capture',{}, ...
    'fullpath',{}, ...
    'rpm',{}, ...
    'bpm',{}, ...
    'rr_accepted',{}, ...
    'hr_accepted',{}, ...
    'status',{}, ...
    'error',{}, ...
    'info',{});

%% 6. PROCESAR TODAS LAS CAPTURAS

fprintf('\nIniciando procesamiento...\n');

for k = 1:N

    archivo = validFiles{k};

    fprintf('\n--------------------------------------------\n');
    fprintf('CAPTURA %d / %d\n',k,N);
    fprintf('%s\n',archivo);
    fprintf('--------------------------------------------\n');

    rpm = NaN;
    bpm = NaN;
    info = struct();
    status = 'OK';
    errorMsg = '';

    try
        % Ejecutar pipeline actual sin generar figuras.
        [rpm,bpm,info] = replotear(archivo,false);

    catch ME
        status = 'ERROR';
        errorMsg = ME.message;

        fprintf(2,'ERROR: %s\n',ME.message);
    end

    % Nombre relativo de la captura.
    relativeName = strrep(archivo, ...
        [captureFolder filesep],'');

    results(end+1) = struct( ...
        'capture',relativeName, ...
        'fullpath',archivo, ...
        'rpm',rpm, ...
        'bpm',bpm, ...
        'rr_accepted',isfinite(rpm), ...
        'hr_accepted',isfinite(bpm), ...
        'status',status, ...
        'error',errorMsg, ...
        'info',info);

    fprintf('RR: %.3f RPM | HR: %.3f BPM\n',rpm,bpm);
end

%% 7. CONSTRUIR TABLA PRINCIPAL

capture = string({results.capture})';
RR_RPM = [results.rpm]';
HR_BPM = [results.bpm]';

RR_accepted = [results.rr_accepted]';
HR_accepted = [results.hr_accepted]';

Status = string({results.status})';
Error_message = string({results.error})';

Tsummary = table(capture,RR_RPM,HR_BPM, ...
    RR_accepted,HR_accepted,Status,Error_message);

%% 8. EXTRAER TODOS LOS CAMPOS DE INFO

diagCapture = strings(0,1);
diagField = strings(0,1);
diagValue = strings(0,1);

for k = 1:numel(results)

    infok = results(k).info;

    if ~isstruct(infok) || isempty(infok)
        continue;
    end

    fields = fieldnames(infok);

    for j = 1:numel(fields)

        field = fields{j};
        value = infok.(field);

        diagCapture(end+1,1) = string(results(k).capture);
        diagField(end+1,1) = string(field);
        diagValue(end+1,1) = valueToText(value);

    end
end

Tdiagnostics = table( ...
    diagCapture,diagField,diagValue, ...
    'VariableNames',{'Capture','Diagnostic','Value'});

%% 9. TABLA DE ARCHIVOS OMITIDOS

Tskipped = table( ...
    string(skippedFiles(:)), ...
    string(skippedReasons(:)), ...
    'VariableNames',{'File','Reason'});

%% 10. EXPORTAR RESULTADOS

xlsxFile = fullfile(outputFolder, ...
    'ECHO_resultados_lote.xlsx');

matFile = fullfile(outputFolder, ...
    'ECHO_resultados_completos.mat');

csvFile = fullfile(outputFolder, ...
    'ECHO_resumen.csv');

writetable(Tsummary,xlsxFile,'Sheet','Summary');

writetable(Tdiagnostics,xlsxFile,'Sheet','Diagnostics');

writetable(Tskipped,xlsxFile,'Sheet','Skipped');

writetable(Tsummary,csvFile);

save(matFile, ...
    'results','Tsummary','Tdiagnostics','Tskipped', ...
    'captureFolder','timestamp','-v7.3');

%% 11. RESUMEN FINAL

fprintf('\n============================================\n');
fprintf('         PROCESAMIENTO COMPLETADO\n');
fprintf('============================================\n');

fprintf('Capturas procesadas: %d\n',N);
fprintf('RR aceptadas: %d / %d\n',sum(RR_accepted),N);
fprintf('HR aceptadas: %d / %d\n',sum(HR_accepted),N);
fprintf('Errores de procesamiento: %d\n', ...
    sum(Status == "ERROR"));

fprintf('\nResultados guardados en:\n%s\n',outputFolder);

fprintf('\nArchivos generados:\n');
fprintf('- ECHO_resultados_lote.xlsx\n');
fprintf('- ECHO_resultados_completos.mat\n');
fprintf('- ECHO_resumen.csv\n');

fprintf('\n============================================\n');

%% FUNCION AUXILIAR

function txt = valueToText(v)

    if isempty(v)
        txt = "<empty>";

    elseif ischar(v)
        txt = string(v);

    elseif isstring(v)
        txt = strjoin(v(:),", ");

    elseif isnumeric(v) || islogical(v)

        if isscalar(v)
            txt = string(sprintf('%.12g',double(v)));
        else
            txt = string(mat2str(v));
        end

    elseif isstruct(v) || iscell(v)

        try
            txt = string(jsonencode(v));
        catch
            txt = "<complex data>";
        end

    else

        try
            txt = string(v);
        catch
            txt = "<unsupported data>";
        end
    end
end
