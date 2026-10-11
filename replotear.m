function [rpm, bpm, info] = replotear(archivo, rango_s, protocolo)
%REPLOTEAR Redibuja el panel de diagnostico de una captura ya guardada.
%
%   Hace lo mismo que main.m pero leyendo de un .mat en vez de capturar:
%   frontend -> filtros -> procesador -> panel de 9 graficas.
%
%   Uso:
%       replotear('captura.mat')                  % igual que main
%       replotear('captura.mat', [15 30])         % solo un tramo
%       replotear('captura.mat', [], 'pautada')   % con el otro paso-alto
%
%   El 'protocolo' debe coincidir con el que usaste al capturar, o los
%   filtros no seran los mismos y el resultado tampoco.

if nargin < 2, rango_s = []; end
if nargin < 3 || isempty(protocolo), protocolo = 'normal'; end

S = load(archivo);
if isfield(S, 'rx_continuous'), rx = S.rx_continuous; else, rx = S.rx; end
cfg = S.cfg;

fprintf('Redibujando %s (%.1f s)...\n', archivo, numel(rx)/cfg.fs);

[iq, Fs] = frontend_dsp(rx, cfg);

switch protocolo
    case 'normal',  fc_hp = 0.08;
    case 'pautada', fc_hp = 0.04;
    otherwise, error('protocolo debe ser ''normal'' o ''pautada''.');
end
fprintf('Protocolo ''%s'': paso-alto en %.2f Hz.\n', protocolo, fc_hp);

dHP = designfilt('highpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', fc_hp, 'SampleRate', Fs);
dLP_rr = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 0.8, 'SampleRate', Fs);

if isempty(rango_s)
    [rpm, bpm, info] = procesar_signos_vitales_fixed(iq, Fs, dHP, dLP_rr);
else
    [rpm, bpm, info] = procesar_signos_vitales_fixed(iq, Fs, dHP, dLP_rr, rango_s);
end

fprintf('\n=============================\n');
fprintf('Respiracion: %.1f rpm\n', rpm);
fprintf('Ritmo cardiaco: %.1f bpm (acuerdo %.0f%%)\n', bpm, 100*info.hr_agreement);
fprintf('=============================\n');

plot_diagnostico_radar(rx, cfg, iq, Fs, rpm, bpm, info);
end
