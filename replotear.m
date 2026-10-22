    function [rpm, bpm, info] = replotear(archivo)
%REPLOTEAR Reprocesa una captura .mat usando el pipeline actual.
%
% Uso:
%   [rpm,bpm,info] = replotear('captura.mat');
%
% El .mat debe contener:
%   rx_continuous (o rx)
%   cfg

S = load(archivo);

% ---------------------------------------------------------
% 1. CARGAR CAPTURA CRUDA
% ---------------------------------------------------------
if isfield(S,'rx_continuous')
    rx = S.rx_continuous;
elseif isfield(S,'rx')
    rx = S.rx;
else
    error('El archivo no contiene rx_continuous ni rx.');
end

if ~isfield(S,'cfg')
    error('El archivo no contiene cfg.');
end

cfg = S.cfg;

if ~isfield(cfg,'fs') || ~isfield(cfg,'tx_tone_hz')
    error('cfg debe contener fs y tx_tone_hz.');
end

fprintf('\n========================================\n');
fprintf('REPROCESANDO: %s\n', archivo);
fprintf('Duracion cruda: %.1f s\n', numel(rx)/cfg.fs);
fprintf('========================================\n');

% ---------------------------------------------------------
% 2. LOW-IF -> BANDA BASE
% ---------------------------------------------------------
fprintf('\nEjecutando front-end DSP...\n');

[iq, Fs] = frontend_dsp(rx, cfg);

% ---------------------------------------------------------
% 3. FILTROS FISIOLOGICOS ACTUALES
% ---------------------------------------------------------
fc_hp = 0.04;

dHP = designfilt('highpassiir', ...
    'FilterOrder',4, ...
    'HalfPowerFrequency',fc_hp, ...
    'SampleRate',Fs);

dLP_rr = designfilt('lowpassiir', ...
    'FilterOrder',4, ...
    'HalfPowerFrequency',0.8, ...
    'SampleRate',Fs);

fprintf('Filtro fisiologico: HP %.2f Hz.\n',fc_hp);

% ---------------------------------------------------------
% 4. PROCESAMIENTO ACTUAL
% ---------------------------------------------------------
[rpm,bpm,info] = procesar_signos_vitales_fixed( ...
    iq, Fs, dHP, dLP_rr);

% ---------------------------------------------------------
% 5. RESULTADOS
% ---------------------------------------------------------
fprintf('\n=============================\n');
fprintf('RESULTADOS\n');

if ~isfinite(rpm)
    if isfield(info,'rr_candidate_rpm') && ...
            isfinite(info.rr_candidate_rpm)

        fprintf(['Respiracion: SIN DETECCION ' ...
            '(candidato %.1f rpm rechazado)\n'], ...
            info.rr_candidate_rpm);
    else
        fprintf('Respiracion: SIN DETECCION\n');
    end
else
    fprintf('Respiracion: %.1f rpm\n',rpm);
end

if ~isfinite(bpm)
    if isfield(info,'hr_candidate_bpm') && ...
            isfinite(info.hr_candidate_bpm)

        fprintf(['Ritmo cardiaco: SIN DETECCION ' ...
            '(candidato %.1f bpm rechazado)\n'], ...
            info.hr_candidate_bpm);
    else
        fprintf('Ritmo cardiaco: SIN DETECCION\n');
    end
else
    fprintf('Ritmo cardiaco: %.1f bpm\n',bpm);
end

if isfield(info,'hr_agreement') && isfinite(info.hr_agreement)
    fprintf('Acuerdo HR entre ventanas: %.0f%%\n', ...
        100*info.hr_agreement);
end

fprintf('=============================\n');

% Regenerar exactamente el mismo diagnostico visual a partir de la captura.
% La figura no interviene en los calculos: solo representa los resultados
% que se acaban de obtener con el pipeline actual.
plot_diagnostico_radar(rx, cfg, iq, Fs, rpm, bpm, info);

end
