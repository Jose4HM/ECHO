function plot_diagnostico_radar(rx_raw, cfg, iq_bb, Fs_bb, rpm, bpm, info)
%PLOT_DIAGNOSTICO_RADAR Panel de 9 graficas de diagnostico de la cadena.
%
%   Estaba como funcion LOCAL dentro de main.m, asi que solo se podia usar
%   desde una captura en vivo. Extraida a archivo propio para poder
%   redibujar capturas guardadas (ver replotear.m).
%
%   Usa las senales que el procesador realmente empleo (info), sin
%   recalcular nada: las graficas y el resultado no pueden discrepar.

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
title(sprintf('Cardiaca tras cancelar respiracion (%.1f bpm)', bpm));

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
    yline(bpm+5, ':k'); yline(bpm-5, ':k');
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
