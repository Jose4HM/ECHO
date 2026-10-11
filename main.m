% =========================================================================
% MASTER SCRIPT: RADAR DE SIGNOS VITALES  (v4)

clear; clc; close all;

% --- PASO 1: CAPTURA ---
try
    [rx_continuous, cfg] = cw();
catch ME
    error(['Fallo en cw.m durante la captura: %s\n' ...
        'Revisa conexion del hardware y configuracion antes de reintentar.'], ...
        ME.message);
end

if ~exist('rx_continuous', 'var') || isempty(rx_continuous)
    error('cw.m no genero rx_continuous.');
end

fprintf('\nEjecutando Front-End DSP...\n');
[iq_frontend, Fs_frontend] = frontend_dsp(rx_continuous, cfg);

if isempty(iq_frontend) || ~isfinite(Fs_frontend) || Fs_frontend <= 0
    error('frontend_dsp devolvio una salida invalida.');
end

% --- PASO 2: FILTROS FISIOLOGICOS ---
% A Fs = 100 Hz estas frecuencias normalizadas quedan en el rango sano
% (0.0016 y 0.016), asi que designfilt produce filtros estables.
%
% EL CORTE DEL PASO-ALTO VA ATADO AL PROTOCOLO DE RESPIRACION:
%
%   'normal'   respiracion espontanea, 12-20 rpm = 0.20-0.33 Hz.
%              Corte en 0.08 Hz, captura de 15-60 s. Es el caso de siempre.
%
%   'pautada'  respiracion guiada con metronomo a 6 rpm = 0.10 Hz, para
%              empujar los armonicos que contaminan la banda cardiaca hasta
%              ordenes altos (el 14o en vez del 3o o 4o), donde su amplitud
%              ya es despreciable. Con corte en 0.08 Hz el filtro se comeria
%              la propia respiracion, asi que baja a 0.04 Hz.
%              PERO un paso-alto de 0.04 Hz tiene ~4 s de transitorio: la
%              captura DEBE ser de 30 s o mas (cfg.capture_seconds en cw.m),
%              o el recorte de bordes se lleva medio registro.
protocolo = 'normal';        % <-- 'normal' | 'pautada'

fc_hp = 0.08;

dur_s = numel(iq_frontend) / Fs_frontend;
if fc_hp < 0.06 && dur_s < 25
    warning(['Protocolo ''%s'' (corte %.2f Hz) con una captura de solo ' ...
        '%.0f s. El transitorio del filtro es de ~%.0f s. Sube ' ...
        'cfg.capture_seconds a 30 en cw.m.'], protocolo, fc_hp, dur_s, ...
        1/(2*pi*fc_hp));
end

fprintf('Protocolo ''%s'': paso-alto en %.2f Hz sobre %.0f s de captura.\n', ...
    protocolo, fc_hp, dur_s);

dHP = designfilt('highpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', fc_hp, 'SampleRate', Fs_frontend);

dLP_rr = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 0.8, 'SampleRate', Fs_frontend);

% --- PASO 3: PROCESAMIENTO ---
fprintf('Datos capturados. Procesando signos vitales...\n');
[rpm, bpm, vital_info] = procesar_signos_vitales_fixed( ...
    iq_frontend, Fs_frontend, dHP, dLP_rr);

% --- PASO 4: RESULTADOS ---
fprintf('\n=============================\n');
fprintf('RESULTADOS FINALES:\n');
if vital_info.rr_fallback
    fprintf('Respiracion: SIN DETECCION\n');
else
    fprintf('Respiracion: %.1f rpm\n', rpm);
end
if vital_info.hr_fallback
    fprintf('Ritmo cardiaco: SIN DETECCION\n');
else
    fprintf('Ritmo cardiaco: %.1f bpm  (acuerdo entre ventanas: %.0f%%)\n', ...
        bpm, 100*vital_info.hr_agreement);
end
fprintf('=============================\n');

if vital_info.rr_fallback
    fprintf(2, ['AVISO: RR = %.1f rpm es un valor GENERICO de fallback, no ' ...
        'una deteccion real. No hubo ningun pico valido en la banda de ' ...
        'respiracion. NO lo presentes como medicion.\n'], rpm);
end
if vital_info.hr_fallback
    fprintf(2, ['AVISO: HR = %.1f bpm es un valor GENERICO de fallback, no ' ...
        'una deteccion real. NO lo presentes como medicion.\n'], bpm);
elseif vital_info.hr_low_confidence
    fprintf(2, ['AVISO: HR = %.1f bpm proviene de un pico real pero de baja ' ...
        'confianza (prominencia %.2f, ratio %.2f, acuerdo %.0f%%). ' ...
        'Tratar con cautela.\n'], bpm, vital_info.hr_prominence, ...
        vital_info.hr_score_ratio, 100*vital_info.hr_agreement);
end

% --- PASO 5: GRAFICAS DE DIAGNOSTICO ---
plot_diagnostico_radar(rx_continuous, cfg, iq_frontend, Fs_frontend, ...
    rpm, bpm, vital_info);


% =========================================================================
function plot_diagnostico_radar(rx_raw, cfg, iq_bb, Fs_bb, rpm, bpm, info)
%PLOT_DIAGNOSTICO_RADAR Resume la cadena usando las senales que el
%   procesador realmente uso (info), sin recalcular nada.

Fs_raw = cfg.fs;
rx_raw = rx_raw(:);
iq_bb  = iq_bb(:);

figure('Name', 'Diagnostico radar signos vitales', ...
    'Units', 'normalized', 'Position', [0.04 0.06 0.92 0.86]);

sgtitle(sprintf(['Diagnostico radar | RR %.1f rpm | HR %.1f bpm | ' ...
    '|c|/R = %.0f | arco %.2f rad'], rpm, bpm, ...
    info.center_offset_ratio, info.arc_span_rad), ...
    'FontWeight', 'bold', 'FontSize', 13);

% ---- 1) Envolvente RX por bloque a lo largo de toda la captura --------
subplot(3, 3, 1);
if isfield(cfg, 'block_peak') && isfield(cfg, 'block_rms')
    n_blk = numel(cfg.block_peak);
    t_blk = (0:n_blk-1) * (cfg.block_len / cfg.fs);
    plot(t_blk, cfg.block_peak, 'b'); hold on;
    plot(t_blk, cfg.block_rms, 'c');
    yline(0.99, ':r');
    idx_inc = cfg.actual_counts ~= cfg.block_len;
    leg = {'pico |bloque|', 'RMS bloque', 'saturacion'};
    if any(idx_inc)
        plot(t_blk(idx_inc), cfg.block_peak(idx_inc), 'rx', ...
            'MarkerSize', 6, 'LineWidth', 1.2);
        leg{end+1} = 'bloque incompleto';
    end
    grid on; xlabel('s'); ylabel('amplitud');
    title(sprintf('Envolvente RX (%d/%d bloques incompletos)', nnz(idx_inc), n_blk));
    legend(leg, 'Location', 'best');
else
    n_raw = min(numel(rx_raw), round(0.01 * Fs_raw));
    t_ms = (0:n_raw-1) / Fs_raw * 1e3;
    plot(t_ms, real(rx_raw(1:n_raw)), 'b'); hold on;
    plot(t_ms, imag(rx_raw(1:n_raw)), 'r');
    grid on; xlabel('ms'); ylabel('amplitud');
    title('RX cruda (cw.m sin instrumentar)');
end

% ---- 2) Constelacion IQ SIN centrar: aqui se ve el arco y su centro ----
subplot(3, 3, 2);
idx_c = round(linspace(1, numel(iq_bb), min(8000, numel(iq_bb))));
plot(real(iq_bb(idx_c)), imag(iq_bb(idx_c)), '.', 'MarkerSize', 3); hold on;
plot(0, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
axis equal; grid on; xlabel('I'); ylabel('Q');
title(sprintf('Constelacion IQ @%.0f Hz (+ = origen)', Fs_bb));

% ---- 3) FFT cruda: debe verse el tono low-IF -------------------------
subplot(3, 3, 3);
plot_fft_complex(rx_raw, Fs_raw, 262144, [-150 150], 1e3);
xline(cfg.tx_tone_hz/1e3, '--r', sprintf('%.0f kHz', cfg.tx_tone_hz/1e3), ...
    'LabelVerticalAlignment', 'bottom');
grid on; xlabel('kHz'); ylabel('dB');
title('FFT RX cruda antes del frontend');

% ---- 4) Fase demodulada en el tiempo ---------------------------------
subplot(3, 3, 4);
t_phi = (0:numel(info.phi)-1).' / info.Fs;
plot(t_phi, info.phi, 'LineWidth', 0.8);
grid on; xlabel('s'); ylabel('rad');
title('Fase desenvuelta (centro del circulo restado)');

% ---- 5) Senal de respiracion en el tiempo ----------------------------
subplot(3, 3, 5);
t_x = (0:numel(info.x_rr)-1).' / info.Fs;
plot(t_x, info.x_rr, 'LineWidth', 1.0);
grid on; xlabel('s'); ylabel('rad');
title(sprintf('Respiracion filtrada (%.1f rpm)', rpm));

% ---- 6) Senal cardiaca en el tiempo ----------------------------------
subplot(3, 3, 6);
t_h = (0:numel(info.x_hr)-1).' / info.Fs;
plot(t_h, info.x_hr, 'LineWidth', 0.8);
grid on; xlabel('s'); ylabel('rad');
title(sprintf('Cardiaca filtrada (%.1f bpm)', bpm));

% ---- 7) FFT respiracion ----------------------------------------------
subplot(3, 3, 7);
plot(info.f_rr_axis, 20*log10(info.P_rr + 1e-12), 'LineWidth', 1.1); hold on;
xline(info.f_rr_est, '--r', sprintf('%.1f rpm', rpm), 'LabelVerticalAlignment', 'bottom');
xline(0.10, ':k', '6 rpm');
xline(0.75, ':k', '45 rpm');
grid on; xlim([0 1.2]); xlabel('Hz'); ylabel('dB');
title('FFT respiracion');

% ---- 8) FFT cardiaca --------------------------------------------------
subplot(3, 3, 8);
plot(info.f_hr_axis, 20*log10(info.P_hr + 1e-12), 'LineWidth', 1.1); hold on;
xline(info.f_hr_est, '--r', sprintf('%.1f bpm', bpm), 'LabelVerticalAlignment', 'bottom');
% Marcar donde caen los armonicos de respiracion: es lo que hay que evitar
for h = 2:6
    fh = h * info.f_rr_est;
    if fh > 0.4 && fh < 3.6
        xline(fh, ':', sprintf('%dxRR', h), 'Color', [0.6 0.6 0.6]);
    end
end
grid on; xlim([0.4 3.6]); xlabel('Hz'); ylabel('dB');
title('FFT cardiaca (gris = armonicos de RR)');

% ---- 9) Consistencia de HR por ventanas -------------------------------
subplot(3, 3, 9);
if ~isempty(info.hr_window_bpm)
    stem(1:numel(info.hr_window_bpm), info.hr_window_bpm, 'filled'); hold on;
    yline(bpm, '--r');
    yline(bpm+5, ':y'); yline(bpm-5, ':y');
    grid on; xlabel('ventana'); ylabel('bpm');
    ylim([max(30, bpm-45) min(200, bpm+45)]);
    title(sprintf('HR por ventana (acuerdo %.0f%%)', 100*info.hr_agreement));
else
    text(0.5, 0.5, 'sin ventanas', 'HorizontalAlignment', 'center');
    axis off;
end
end

% =========================================================================
function plot_fft_complex(x, Fs, nfft_max, xlim_values, unit_scale)
x = x(:);
N = min(numel(x), nfft_max);
Nfft = 2^nextpow2(N);
w = hamming(N);
X = fftshift(fft(x(1:N) .* w, Nfft));
f = (-Nfft/2:Nfft/2-1) * Fs / Nfft;
plot(f / unit_scale, 20*log10(abs(X) + 1e-12), 'LineWidth', 1.0);
xlim(xlim_values);
end
