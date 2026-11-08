function calibrar_antenas_live()
%CALIBRAR_ANTENAS_LIVE Monitor RF en vivo para orientar TX/RX del radar CW.
%
% Uso:
%   calibrar_antenas_live
%
% Transmite continuamente un tono Low-IF de 100 kHz por TX1 y recibe por
% RX2. Muestra espectro, trayectoria IQ, nivel del tono, tono/piso y fase.
%
% IMPORTANTE:
%   - Cerrar la ventana detiene la calibracion.
%   - No maximizar solamente el nivel recibido: con TX/RX RHCP puede dominar
%     el acoplamiento directo entre antenas.
%   - Para orientar al pecho, buscar una recepcion estable, sin saturacion,
%     y una variacion IQ/fase visible cuando el blanco se mueve o respira.

clc;
fprintf('==============================================\n');
fprintf(' CALIBRACION EN VIVO DE ANTENAS - RADAR CW\n');
fprintf('==============================================\n');

%% CONFIGURACION
cfg = struct();
cfg.serial       = '8c14dec0ebc947b88424691d0dbf744f';
cfg.fc           = 5.94e9;
cfg.fs           = 1e6;
cfg.bandwidth    = 1e6;
cfg.tx_tone_hz   = 100e3;
cfg.tx_amplitude = 0.25;
cfg.tx_channel   = 'TX1';
cfg.rx_channel   = 'RX2';
cfg.tx_gain      = 65;
cfg.rx_gain      = 0;
cfg.rx_agc       = 'manual';

cfg.block_len     = 30720;
cfg.num_buffers   = 32;
cfg.buffer_size   = 30720;
cfg.num_transfers = 8;
cfg.timeout_ms    = 5000;

cfg.display_span_hz = 20e3;
cfg.update_every    = 4;
cfg.history_s       = 20;
cfg.motion_window_s = 2.0;

assert(cfg.tx_amplitude > 0 && cfg.tx_amplitude <= 0.5, ...
    'tx_amplitude debe estar entre 0 y 0.5.');

samples_per_period = cfg.fs / cfg.tx_tone_hz;
assert(abs(samples_per_period - round(samples_per_period)) < 1e-12, ...
    'fs/tx_tone_hz debe ser entero.');
assert(mod(cfg.block_len, round(samples_per_period)) == 0, ...
    'block_len debe contener un numero entero de periodos del tono.');

if exist('bladeRF', 'class') ~= 8 && exist('bladeRF', 'file') ~= 2
    error('No se encontro bladeRF.m en el path de MATLAB.');
end

%% ABRIR Y CONFIGURAR BLADERF
fprintf('Abriendo bladeRF...\n');
b = bladeRF(['*:serial=' cfg.serial]);
cleanup_obj = onCleanup(@() cleanup_bladerf_local(b, cfg));

b.tx.channel    = cfg.tx_channel;
b.rx.channel    = cfg.rx_channel;
b.tx.frequency  = cfg.fc;
b.rx.frequency  = cfg.fc;
b.tx.samplerate = cfg.fs;
b.rx.samplerate = cfg.fs;
b.tx.bandwidth  = cfg.bandwidth;
b.rx.bandwidth  = cfg.bandwidth;

modo_puesto = '';
for nombre = {'agc', 'gain_mode'}
    campo = nombre{1};
    try
        b.rx.(campo) = cfg.rx_agc;
        modo_puesto = campo;
        fprintf('RX: %s = %s\n', campo, string(b.rx.(campo)));
        break;
    catch
    end
end
if isempty(modo_puesto)
    warning(['No se pudo fijar explicitamente el modo manual de RX. ' ...
        'Se continua porque este binding puede ignorar esa propiedad.']);
end

b.tx.gain = cfg.tx_gain;
pause(0.05);
try
    tx_gain_real = b.tx.gain;
catch
    tx_gain_real = NaN;
end

try
    b.rx.gain = cfg.rx_gain;
    pause(0.05);
    rx_gain_real = b.rx.gain;
catch
    rx_gain_real = NaN;
end

fprintf('TX gain: solicitado %.1f dB | reportado %.1f dB\n', ...
    cfg.tx_gain, tx_gain_real);
fprintf('RX gain: solicitado %.1f dB | reportado %.1f dB\n', ...
    cfg.rx_gain, rx_gain_real);

b.tx.config.num_buffers   = cfg.num_buffers;
b.tx.config.buffer_size   = cfg.buffer_size;
b.tx.config.num_transfers = cfg.num_transfers;
b.tx.config.timeout_ms    = cfg.timeout_ms;

b.rx.config.num_buffers   = cfg.num_buffers;
b.rx.config.buffer_size   = cfg.buffer_size;
b.rx.config.num_transfers = cfg.num_transfers;
b.rx.config.timeout_ms    = cfg.timeout_ms;

%% TONO TX
n = (0:cfg.block_len-1).';
tone = exp(1j*2*pi*cfg.tx_tone_hz*n/cfg.fs);
tx_block = cfg.tx_amplitude * tone;
tx_zeros = complex(zeros(cfg.block_len,1));

fprintf('Iniciando RX %s y TX %s...\n', cfg.rx_channel, cfg.tx_channel);
b.rx.start();
b.tx.start();

fprintf('Inicializando streaming...\n');
for k = 1:20
    b.transmit(tx_block, cfg.timeout_ms, 0, false, false);
    [~, ~, ~, ~] = b.receive(cfg.block_len);
end

%% FIGURA
fig = figure( ...
    'Name', 'Calibracion de antenas - radar CW', ...
    'NumberTitle', 'off', ...
    'Units', 'normalized', ...
    'Position', [0.05 0.08 0.90 0.82]);

tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axSpec = nexttile(tl, 1);
hSpec = plot(axSpec, nan, nan, 'LineWidth', 1.0);
hold(axSpec, 'on');
xline(axSpec, cfg.tx_tone_hz/1e3, '--', 'TX nominal');
grid(axSpec, 'on');
xlabel(axSpec, 'Frecuencia (kHz)');
ylabel(axSpec, 'Amplitud (dBFS aprox.)');
title(axSpec, 'Espectro RX alrededor del tono');
xlim(axSpec, [(cfg.tx_tone_hz-cfg.display_span_hz)/1e3, ...
              (cfg.tx_tone_hz+cfg.display_span_hz)/1e3]);
ylim(axSpec, [-120 5]);

axIQ = nexttile(tl, 2);
hIQ = plot(axIQ, nan, nan, '.-', 'MarkerSize', 6);
hold(axIQ, 'on');
hIQLast = plot(axIQ, nan, nan, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
grid(axIQ, 'on');
axis(axIQ, 'equal');
xlabel(axIQ, 'I coherente');
ylabel(axIQ, 'Q coherente');
title(axIQ, 'IQ del tono recibido - historial');

axLevel = nexttile(tl, 3);
yyaxis(axLevel, 'left');
hToneHist = plot(axLevel, nan, nan, 'LineWidth', 1.1);
ylabel(axLevel, 'Nivel del tono (dBFS aprox.)');
ylim(axLevel, [-100 5]);
yyaxis(axLevel, 'right');
hSNRHist = plot(axLevel, nan, nan, 'LineWidth', 1.1);
ylabel(axLevel, 'Tono / piso (dB)');
ylim(axLevel, [0 100]);
grid(axLevel, 'on');
xlabel(axLevel, 'Tiempo (s)');
title(axLevel, 'Nivel recibido y separacion respecto al piso');

axPhase = nexttile(tl, 4);
hPhase = plot(axPhase, nan, nan, 'LineWidth', 1.1);
grid(axPhase, 'on');
xlabel(axPhase, 'Tiempo (s)');
ylabel(axPhase, 'Fase relativa (rad)');
title(axPhase, 'Fase del tono recibido');

sgtitle(tl, 'Esperando datos...', 'FontWeight', 'bold');

%% VARIABLES DE MONITOREO
block_dt = cfg.block_len / cfg.fs;
max_hist = max(20, ceil(cfg.history_s / block_dt));

t_hist      = zeros(0,1);
z_hist      = complex(zeros(0,1));
tone_hist   = zeros(0,1);
snr_hist    = zeros(0,1);

Nfft = 2^nextpow2(cfg.block_len);
w = hann(cfg.block_len, 'periodic');
w_sum = sum(w);
f_ax = (-Nfft/2:Nfft/2-1).' * cfg.fs/Nfft;

plot_band = abs(f_ax - cfg.tx_tone_hz) <= cfg.display_span_hz;
tone_band = abs(f_ax - cfg.tx_tone_hz) <= 500;
noise_band = abs(f_ax - cfg.tx_tone_hz) >= 5e3 & ...
             abs(f_ax - cfg.tx_tone_hz) <= cfg.display_span_hz;

total_overruns = 0;
block_counter = 0;
sample_counter = 0;

fprintf('\nCalibracion activa.\n');
fprintf('Cierra la ventana para terminar.\n');
fprintf(['No busques solo el maximo nivel: una posicion util debe evitar ' ...
         'saturacion y mostrar respuesta IQ/fase al movimiento del blanco.\n\n']);

%% LOOP EN VIVO
try
    while isgraphics(fig)
        block_counter = block_counter + 1;

        b.transmit(tx_block, cfg.timeout_ms, 0, false, false);
        [rx, ~, actual_count, overrun] = b.receive(cfg.block_len);

        total_overruns = total_overruns + double(overrun ~= 0);

        if actual_count <= 0
            drawnow;
            continue;
        end

        rx = rx(1:actual_count);
        nloc = (0:actual_count-1).';

        % Coeficiente complejo coherente en el tono nominal.
        lo = exp(-1j*2*pi*cfg.tx_tone_hz*nloc/cfg.fs);
        z = mean(rx .* lo);

        % Pico por componente ADC; |I+jQ| no se usa como criterio de clipping.
        peak_adc = max([abs(real(rx)); abs(imag(rx))]);

        if actual_count == cfg.block_len
            X = fftshift(fft(rx .* w, Nfft));
            Pdb = 20*log10(abs(X)/(w_sum + eps) + 1e-12);

            tone_db = max(Pdb(tone_band));
            noise_db = median(Pdb(noise_band));
            tone_to_floor_db = tone_db - noise_db;
        else
            tone_db = 20*log10(abs(z) + 1e-12);
            tone_to_floor_db = NaN;
            Pdb = [];
        end

        sample_counter = sample_counter + actual_count;
        t_now = sample_counter / cfg.fs;

        t_hist(end+1,1)    = t_now; %#ok<AGROW>
        z_hist(end+1,1)    = z; %#ok<AGROW>
        tone_hist(end+1,1) = tone_db; %#ok<AGROW>
        snr_hist(end+1,1)  = tone_to_floor_db; %#ok<AGROW>

        if numel(t_hist) > max_hist
            keep = (numel(t_hist)-max_hist+1):numel(t_hist);
            t_hist    = t_hist(keep);
            z_hist    = z_hist(keep);
            tone_hist = tone_hist(keep);
            snr_hist  = snr_hist(keep);
        end

        if mod(block_counter, cfg.update_every) ~= 0
            continue;
        end

        if ~isempty(Pdb)
            set(hSpec, 'XData', f_ax(plot_band)/1e3, ...
                       'YData', Pdb(plot_band));
        end

        set(hIQ, 'XData', real(z_hist), 'YData', imag(z_hist));
        set(hIQLast, 'XData', real(z_hist(end)), 'YData', imag(z_hist(end)));
        auto_limits_iq(axIQ, z_hist);

        t_rel = t_hist - t_hist(end);
        yyaxis(axLevel, 'left');
        set(hToneHist, 'XData', t_rel, 'YData', tone_hist);
        yyaxis(axLevel, 'right');
        set(hSNRHist, 'XData', t_rel, 'YData', snr_hist);
        xlim(axLevel, [-cfg.history_s 0]);

        ph = unwrap(angle(z_hist));
        ph = ph - ph(1);
        set(hPhase, 'XData', t_rel, 'YData', ph);
        xlim(axPhase, [-cfg.history_s 0]);

        n_motion = max(3, round(cfg.motion_window_s / block_dt));
        idxm = max(1, numel(ph)-n_motion+1):numel(ph);
        ph_recent = ph(idxm);
        if numel(ph_recent) >= 3
            ph_dyn = detrend(ph_recent, 1);
            phase_pp = max(ph_dyn) - min(ph_dyn);
        else
            phase_pp = NaN;
        end

        if peak_adc >= 0.98
            estado = 'SATURADO - separa antenas o baja nivel TX';
        elseif peak_adc >= 0.90
            estado = 'CERCA DE SATURACION';
        elseif isfinite(tone_to_floor_db) && tone_to_floor_db < 10
            estado = 'TONO DEBIL / PISO ALTO';
        elseif isfinite(tone_to_floor_db) && tone_to_floor_db >= 20
            estado = 'RF ESTABLE';
        else
            estado = 'RF UTIL - revisar graficas';
        end

        sgtitle(tl, sprintf([ ...
            '%s | tono %.1f dBFS | tono/piso %.1f dB | pico ADC %.3f | ' ...
            'dinamica fase(%.0fs) %.3f rad | overruns %d'], ...
            estado, tone_db, tone_to_floor_db, peak_adc, ...
            cfg.motion_window_s, phase_pp, total_overruns), ...
            'FontWeight', 'bold');

        drawnow limitrate;
    end

catch ME
    fprintf(2, '\nError durante calibracion: %s\n', ME.message);
    rethrow(ME);
end

try
    b.transmit(tx_zeros, cfg.timeout_ms, 0, false, true);
catch
end

clear cleanup_obj;
fprintf('\nCalibracion finalizada.\n');
end

% =========================================================================
function auto_limits_iq(ax, z)
if isempty(z)
    return;
end
xr = real(z);
yi = imag(z);
xc = median(xr);
yc = median(yi);
span = max([max(xr)-min(xr), max(yi)-min(yi), 1e-4]);
half = 0.65 * span;
xlim(ax, [xc-half, xc+half]);
ylim(ax, [yc-half, yc+half]);
end

% =========================================================================
function cleanup_bladerf_local(b, cfg)
try
    z = complex(zeros(cfg.block_len,1));
    b.transmit(z, cfg.timeout_ms, 0, false, true);
catch
end
try
    b.tx.stop();
catch
end
pause(0.05);
try
    b.rx.stop();
catch
end
end
