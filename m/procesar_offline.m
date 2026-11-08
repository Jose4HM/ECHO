function [rpm, bpm, info] = procesar_offline(archivo, rango_s)
%PROCESAR_OFFLINE Alias compatible de REPLOTEAR.
% Conserva el nombre antiguo, pero delega en la unica ruta validada. No se
% mantiene una segunda configuracion de filtros que pueda cambiar resultados.

if nargin < 1
    error('Uso: procesar_offline(''captura_01.mat'')');
end
if nargin >= 2 && ~isempty(rango_s)
    error(['El procesamiento por tramo esta deshabilitado porque usaba una ' ...
        'ruta distinta al pipeline validado. Reprocese la captura completa.']);
end

[rpm, bpm, info] = replotear(archivo);
end
