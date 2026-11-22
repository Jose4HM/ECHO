%% ECHO - Recalculo auditable de metricas y estudio de ablacion
% Requiere en la misma carpeta:
%  1) resultados_ablacion_ECHO_revisado.xlsx (hoja Resultados)
%  2) referencias_tabla10_ECHO.csv (70 referencias extraidas de la Tabla 10)
%
% NO procesa los archivos MAT ni cambia los algoritmos. Recalcula las
% metricas a partir de las salidas detalladas YA obtenidas en cada variante.
% No excluye valores finitos por tener errores altos.
% Calcula por separado HR y RR y sus subconjuntos comunes de aceptacion.
% Compatible con MATLAB reciente (readtable, writetable, unique, ismember).

clear; clc;

archivo_ablacion = 'resultados_ablacion_ECHO_revisado.xlsx';
archivo_ref = 'referencias_tabla10_ECHO.csv';

assert(isfile(archivo_ablacion), 'No existe %s', archivo_ablacion);
assert(isfile(archivo_ref), 'No existe %s', archivo_ref);

A = readtable(archivo_ablacion, 'Sheet', 'Resultados', ...
    'TextType', 'string', 'VariableNamingRule', 'preserve');
R = readtable(archivo_ref, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');

requeridasA = {'Archivo','Variante','HR_bpm','RR_rpm'};
requeridasR = {'Archivo','HR_ref_bpm','RR_ref_rpm', ...
    'HR_ECHO_tabla10_bpm','RR_ECHO_tabla10_rpm'};
assert(all(ismember(requeridasA, A.Properties.VariableNames)), ...
    'Faltan columnas en la hoja Resultados del archivo de ablacion.');
assert(all(ismember(requeridasR, R.Properties.VariableNames)), ...
    'Faltan columnas en el CSV de referencias.');

R.Archivo = string(R.Archivo);
A.Archivo = string(A.Archivo);
A.Variante = string(A.Variante);
assert(height(R)==70 && numel(unique(R.Archivo))==70, ...
    'El CSV de Tabla 10 debe contener exactamente 70 identificadores unicos.');

% Orden fijo de capturas y variantes (incluye control del procesamiento
% original y copia de ablacion completa para comprobar equivalencia).
ids = R.Archivo;
variantes = unique(A.Variante, 'stable');
nC = numel(ids);
nV = numel(variantes);

HRref = a_num(R.HR_ref_bpm);
RRref = a_num(R.RR_ref_rpm);
HRtabla = a_num(R.HR_ECHO_tabla10_bpm);
RRtabla = a_num(R.RR_ECHO_tabla10_rpm);
assert(all(isfinite(HRref)) && all(isfinite(RRref)), ...
    'La Tabla 10 debe contener las 70 referencias HR y RR numericas.');

HR = nan(nC,nV);
RR = nan(nC,nV);
for v=1:nV
    filas = find(A.Variante==variantes(v));
    assert(numel(filas)==nC, ...
        'Variante %s: se esperaban 70 filas, encontradas %d.', ...
        variantes(v), numel(filas));
    archivosV = A.Archivo(filas);
    assert(numel(unique(archivosV))==nC, ...
        'La variante %s contiene capturas duplicadas.', variantes(v));
    [encontrado, posicion] = ismember(ids, archivosV);
    assert(all(encontrado), 'Faltan capturas de Tabla 10 en %s.', variantes(v));
    HR(:,v) = a_num(A.HR_bpm(filas(posicion)));
    RR(:,v) = a_num(A.RR_rpm(filas(posicion)));
end

% -------------------------------------------------------------------------
% 1) Auditar consistencia entre Tabla 10 (redondeada a un decimal) y el
%    resultado sin redondear del control original. No modificar cifras.
% -------------------------------------------------------------------------
idxOriginal = find(variantes=="Original_control",1);
assert(~isempty(idxOriginal), 'Falta Original_control en el archivo de ablacion.');
fprintf('\n=========== COMPROBACION TABLA 10 ===========\n');
comprobar_tabla(ids,'HR',HRtabla,HR(:,idxOriginal));
comprobar_tabla(ids,'RR',RRtabla,RR(:,idxOriginal));

% -------------------------------------------------------------------------
% 2) Metricas de cada variante sobre su subconjunto propio de salidas finitas
% -------------------------------------------------------------------------
Capturas = repmat(nC,nV,1);
Variante = variantes;
HR_aceptadas = zeros(nV,1);
HR_pares = zeros(nV,1);
HR_cobertura_pct = zeros(nV,1);
HR_MAE = nan(nV,1); HR_RMSE = nan(nV,1);
RR_aceptadas = zeros(nV,1);
RR_pares = zeros(nV,1);
RR_cobertura_pct = zeros(nV,1);
RR_MAE = nan(nV,1); RR_RMSE = nan(nV,1);

for v=1:nV
    [HR_aceptadas(v),HR_pares(v),HR_cobertura_pct(v), ...
        HR_MAE(v),HR_RMSE(v)] = evaluar(HR(:,v),HRref);
    [RR_aceptadas(v),RR_pares(v),RR_cobertura_pct(v), ...
        RR_MAE(v),RR_RMSE(v)] = evaluar(RR(:,v),RRref);
end

Resumen = table(Variante,Capturas,HR_aceptadas,HR_pares, ...
    HR_cobertura_pct,HR_MAE,HR_RMSE,RR_aceptadas,RR_pares, ...
    RR_cobertura_pct,RR_MAE,RR_RMSE);

% -------------------------------------------------------------------------
% 3) Subconjunto COMUN: interseccion de aceptadas en TODAS las variantes.
%    Se determina independientemente para HR y RR.
% -------------------------------------------------------------------------
comunHR = all(isfinite(HR),2) & isfinite(HRref);
comunRR = all(isfinite(RR),2) & isfinite(RRref);
nComunHR = nnz(comunHR);
nComunRR = nnz(comunRR);
assert(nComunHR>0 && nComunRR>0, ...
    'No existen capturas comunes finitas en HR o RR.');

HR_MAE_comun = nan(nV,1); HR_RMSE_comun = nan(nV,1);
RR_MAE_comun = nan(nV,1); RR_RMSE_comun = nan(nV,1);
for v=1:nV
    e = HR(comunHR,v)-HRref(comunHR);
    HR_MAE_comun(v)=mean(abs(e));
    HR_RMSE_comun(v)=sqrt(mean(e.^2));
    e = RR(comunRR,v)-RRref(comunRR);
    RR_MAE_comun(v)=mean(abs(e));
    RR_RMSE_comun(v)=sqrt(mean(e.^2));
end

Comun = table(Variante, repmat(nComunHR,nV,1), ...
    HR_MAE_comun,HR_RMSE_comun, repmat(nComunRR,nV,1), ...
    RR_MAE_comun,RR_RMSE_comun, ...
    'VariableNames', {'Variante','HR_n_comun','HR_MAE_comun', ...
    'HR_RMSE_comun','RR_n_comun','RR_MAE_comun','RR_RMSE_comun'});

% -------------------------------------------------------------------------
% 4) Exportar datos auditables y figuras HR/RR del original (sin NaN)
% -------------------------------------------------------------------------
Detalle = table(ids,HRref,RRref,HRtabla,RRtabla, ...
    comunHR,comunRR, ...
    'VariableNames', {'Archivo','HR_referencia','RR_referencia', ...
    'HR_tabla10_redondeado','RR_tabla10_redondeado', ...
    'Pertenece_comun_HR','Pertenece_comun_RR'});
for v=1:nV
    sufijo = matlab.lang.makeValidName(char(variantes(v)));
    Detalle.(['HR_' sufijo]) = HR(:,v);
    Detalle.(['RR_' sufijo]) = RR(:,v);
end

writetable(Resumen,'ECHO_ablacion_metricas_corregidas.csv');
writetable(Comun,'ECHO_ablacion_comparacion_comun.csv');
writetable(Detalle,'ECHO_ablacion_detalle_70_capturas.csv');

fprintf('\n=========== CADA VARIANTE: TODAS LAS ACEPTADAS ===========\n');
disp(Resumen);
fprintf('\n=========== MISMO SUBCONJUNTO PARA TODAS ===========\n');
fprintf('Capturas comunes: HR=%d de 70; RR=%d de 70\n',nComunHR,nComunRR);
disp(Comun);
fprintf('\nArchivos generados:\n');
fprintf(' ECHO_ablacion_metricas_corregidas.csv\n');
fprintf(' ECHO_ablacion_comparacion_comun.csv\n');
fprintf(' ECHO_ablacion_detalle_70_capturas.csv\n');

%% Funciones auxiliares
function x=a_num(v)
    if isnumeric(v)
        x=double(v);
    else
        x=str2double(string(v));
    end
    x=x(:);
end
function [aceptadas,pares,cobertura,mae,rmse]=evaluar(est,ref)
    aceptadas=nnz(isfinite(est));
    validos=isfinite(est)&isfinite(ref);
    pares=nnz(validos);
    cobertura=100*aceptadas/numel(est);
    mae=NaN;rmse=NaN;
    if pares>0
        err=est(validos)-ref(validos);
        mae=mean(abs(err));
        rmse=sqrt(mean(err.^2));
    end
end
function comprobar_tabla(ids,senal,tbl,calc)
    finTbl=isfinite(tbl); finCalc=isfinite(calc);
    fprintf('%s: numericos tabla=%d, numericos archivo=%d\n', ...
        senal,nnz(finTbl),nnz(finCalc));
    if any(xor(finTbl,finCalc))
        fprintf(' ADVERTENCIA: desacuerdo en aceptacion/rechazo:\n');
        disp(ids(xor(finTbl,finCalc)));
    end
    idx=finTbl & finCalc & (abs(tbl-calc)>0.051);
    if any(idx)
        fprintf(' Diferencias mayores a 0.051 respecto al redondeo a 0.1:\n');
        disp(table(ids(idx),tbl(idx),calc(idx), ...
            'VariableNames',{'Captura','Tabla10','Procesado'}));
    else
        fprintf(' Valores consistentes con redondeo a 0.1.\n');
    end
end
