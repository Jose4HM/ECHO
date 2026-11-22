function [rx_continuous, cfg] = cw()
% =========================================================================
% ECHO - Experimental Continuous-Wave Human Observation
% Copyright (C) 2026 Jose A. Hancco Mamani
% SPDX-License-Identifier: GPL-3.0-only
%
% This program is free software: you can redistribute it
% and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation,
% version 3 of the License.
% =========================================================================

% -------------------------------------------------------------------------
% CW - Continuous-Wave I/Q Data Acquisition
%
% Description:
%   Configures the Nuand bladeRF 2.0 micro xA4 SDR and performs
%   continuous complex I/Q data acquisition for the ECHO platform.
%
%   The function generates a complex Low-IF transmission tone,
%   configures the RF hardware, and records received I/Q samples.
%   It also performs preliminary diagnostics of receiver levels,
%   streaming continuity, incomplete blocks, and signal amplitude.
%
% Syntax:
%   [rx_continuous, cfg] = cw()
%
% Outputs:
%   rx_continuous - Acquired complex I/Q samples.
%   cfg           - Structure containing hardware configuration,
%                   acquisition parameters, and quality diagnostics.
%
% Acquisition stages:
%   1. Hardware configuration and initialization.
%   2. Zero-input TX baseline measurement.
%   3. TX settling period.
%   4. Optional TX amplitude autocalibration.
%   5. Continuous I/Q acquisition.
%   6. Acquisition integrity and signal-level assessment.
%
% Requirements:
%   MATLAB, bladeRF MATLAB interface, and a bladeRF 2.0 micro xA4 SDR.
%
% Note:
%   This function acquires raw complex I/Q samples. Physiological
%   signal processing is performed by subsequent ECHO modules.

cfg = struct();              % Initialize configuration structure for bladeRF parameters
%cfg.serial = '8c14dec0ebc947b88424691d0dbf744f';
cfg.serial = 'YOUR_BLADERF_SERIAL'; % Replace with the serial number of your SDR.

% RF parameters: carrier frequency, sampling rate, and bandwidth.
cfg.fc = 5.94e9;             % RF carrier frequency (Hz).
cfg.fs = 1e6;                % Sampling rate (samples/s).
cfg.bandwidth = 1e6;         % RF bandwidth (Hz).
cfg.tx_tone_hz = 100e3;      % Low-IF transmission tone (Hz).
cfg.tx_amplitude = 0.25;     % Digital TX signal amplitude.

% Hardware configuration: TX/RX channels, gain, and AGC mode.
cfg.tx_channel = 'TX1';      % Transmit channel.
cfg.rx_channel = 'RX2';     % Receive channel.
cfg.tx_gain = 65;            % TX hardware gain (dB).
cfg.rx_gain = 60;            % RX hardware gain (dB).
cfg.rx_agc = 'manual';       % Receiver gain control mode.

% Receiver gain verification settings.
cfg.strict_gain = false;    % Abort acquisition if the applied gain differs excessively.
cfg.gain_tol_db = 3;        % Maximum allowed gain deviation (dB).

cfg.baseline_blocks = 10;
% TX settling period before the main acquisition.
% Initial samples are discarded to allow the transmit chain
% and streaming pipeline to reach a steady operating condition.
cfg.tx_settle_seconds = 20;


% Optional TX digital amplitude autocalibration.
% Disabled by default.
cfg.autocal_enable = false;
cfg.autocal_target_peak = 0.45; % Target received signal peak.
cfg.autocal_max_iter = 6;       % Maximum calibration iterations.
cfg.autocal_blocks = 5;         % Blocks used per calibration measurement.

% Acquisition duration and streaming configuration.
cfg.block_len = 10000;       % Samples per acquisition block.
cfg.capture_seconds = 60;    % Requested acquisition duration (s).

cfg.num_blocks = ceil(cfg.capture_seconds * cfg.fs / cfg.block_len);
cfg.capture_seconds = (cfg.block_len * cfg.num_blocks) / cfg.fs;
% USB buffer and capture duration settings designed to ensure continuous streaming
cfg.num_buffers = 32;
cfg.buffer_size = 30720;
cfg.num_transfers = 8;
cfg.timeout_ms = 5000;

assert(cfg.tx_amplitude > 0 && cfg.tx_amplitude <= 0.5, ...
    'cfg.tx_amplitude must be between 0 and 0.5.');

samples_per_period = cfg.fs / cfg.tx_tone_hz;
assert(abs(samples_per_period - round(samples_per_period)) < 1e-12, ...
    'fs/tx_tone_hz must be an integer to repeat the tone without a phase shift.');
samples_per_period = round(samples_per_period);
assert(mod(cfg.block_len, samples_per_period) == 0, ...
    'cfg.block_len must contain an integer number of periods of the tone.');

if exist('bladeRF', 'class') ~= 8 && exist('bladeRF', 'file') ~= 2
    error('bladeRF.m was not found in the path.');
end

fprintf('Opening bladeRF...\n');
b = bladeRF(['*:serial=' cfg.serial]);
cleanup_obj = onCleanup(@() cleanup_bladerf(b, cfg));

cfg = configure_bladerf(b, cfg);

total_samples = cfg.block_len * cfg.num_blocks;
rx_continuous = complex(zeros(total_samples, 1));
actual_counts = zeros(cfg.num_blocks, 1);
block_peak = zeros(cfg.num_blocks, 1);
block_rms  = zeros(cfg.num_blocks, 1);

% Generate the baseband Low-IF complex tone and zero-array for the transmitter
n = (0:cfg.block_len-1).';
tone = exp(1j*2*pi*cfg.tx_tone_hz*n/cfg.fs);
tx_block = cfg.tx_amplitude * tone;
tx_zeros = complex(zeros(cfg.block_len, 1));

% =========================================================================
% 1. HARDWARE INITIALIZATION AND WARM-UP
% Enable RX/TX, apply the receiver gain, and flush USB buffers
% before starting the main acquisition.
% =========================================================================
fprintf('Waiting for the initial lock of synthesizers to %.3f GHz...\n', cfg.fc/1e9);
pause(1.0);

fprintf('Starting RX %s and TX %s...\n', cfg.rx_channel, cfg.tx_channel);

% 1. Enable receiver
b.rx.start();

% 2. Apply gain with RX enabled
b.rx.gain = cfg.rx_gain;
pause(0.1);

% 3. Read effective gain
cfg.rx_gain_actual = double(b.rx.gain);

fprintf('\n--- RX GAIN VERIFICATION ---\n');
fprintf('Requested: %.1f dB\n', cfg.rx_gain);
fprintf('Applied:   %.1f dB\n', cfg.rx_gain_actual);

% 4. Verify the result
verificar_ganancia_rx(cfg);

% 5. Enable the transmitter.
b.tx.start();

fprintf('Clearing USB Buffers...\n');
for k = 1:30
    b.transmit(tx_zeros, cfg.timeout_ms, 0, false, false);
    [~, ~, ~, ov_warmup] = b.receive(cfg.block_len);
    if ov_warmup
        warning('Overrun during warm-up in an iteration %d', k);
    end
end

% =========================================================================
% 2. ZERO-INPUT TX BASELINE
% Measure the receiver noise floor and external interference
% while the TX stream is driven with zero-valued samples.
% =========================================================================
fprintf('\n--- ZERO-INPUT TX BASELINE ---\n');
[pk_off, rms_off, pk_comp_off] = medir_nivel(b, cfg, tx_zeros, cfg.baseline_blocks);
fprintf('Peak |IQ| without TX: %.4f   I/Q component peak: %.4f\n', pk_off, pk_comp_off);
fprintf('RMS without TX: %.5f (%.1f dBFS)\n', rms_off, 20*log10(rms_off + 1e-12));

% =========================================================================
% 3. TX SETTLING PERIOD AND TX-ON LEVEL ASSESSMENT
% Transmit the Low-IF tone and discard the initial received blocks.
% Then measure the received signal level and compare it with
% the zero-input TX baseline.
% =========================================================================
fprintf('\n--- TX STABILIZATION ON (%.1f s) ---\n', cfg.tx_settle_seconds);
n_settle = ceil(cfg.tx_settle_seconds * cfg.fs / cfg.block_len);
settle_overruns = 0;
settle_incomplete = 0;

% Turn on TX and deliberately discard initial data while the RF/LO phase and TX tail stabilize
for k = 1:n_settle
    b.transmit(tx_block, cfg.timeout_ms, 0, false, false);
    [~, ~, nc_settle, ov_settle] = b.receive(cfg.block_len);
    settle_overruns = settle_overruns + double(ov_settle ~= 0);
    settle_incomplete = settle_incomplete + double(nc_settle ~= cfg.block_len);
end
cfg.tx_settle_actual_seconds = n_settle * cfg.block_len / cfg.fs;
cfg.tx_settle_overruns = settle_overruns;
cfg.tx_settle_incomplete_blocks = settle_incomplete;
fprintf('Stabilization complete: %.1f s, overruns %d, incomplete blocks %d.\n', ...
    cfg.tx_settle_actual_seconds, settle_overruns, settle_incomplete);
if settle_overruns > 0 || settle_incomplete > 0
    warning(['There were issues during the TX stabilization. They are not part ' ...
        'of the usable capture, but it is a good idea to review the stream.']);
end

[pk_on, rms_on, pk_comp_on] = medir_nivel(b, cfg, tx_block, cfg.baseline_blocks);
fprintf('Peak |IQ| with TX: %.4f   I/Q component peak: %.4f\n', pk_on, pk_comp_on);
fprintf('RMS with TX: %.5f (%.1f dBFS)\n', rms_on, 20*log10(rms_on + 1e-12));
razon_db = 20*log10((rms_on + 1e-12) / (rms_off + 1e-12));
fprintf('RMS Increase When TX Is Enabled: %.1f dB\n', razon_db);

cfg.baseline_peak_off = pk_off;
cfg.baseline_rms_off  = rms_off;
cfg.baseline_peak_on  = pk_on;
cfg.baseline_peak_component_off = pk_comp_off;
cfg.baseline_peak_component_on  = pk_comp_on;
cfg.baseline_rms_on   = rms_on;
cfg.leak_over_floor_db = razon_db;

if pk_comp_off >= 0.99
    error(['The receiver is overloaded with TX digitally turned off (peak I/Q %.3f). ' ...
        'Check RX gain and external interference around %.3f GHz before capturing.'], ...
        pk_comp_off, cfg.fc/1e9);
end

if razon_db < 6
    warning(['The RMS increases by only %.1f dB when TX is activated. Verify that the ' ...
        'Low-IF tone is present and that the TX and antennas are functioning properly.'], razon_db);
end

% =========================================================================
% 4. OPTIONAL TX AMPLITUDE AUTOCALIBRATION
% Verify the receiver response to changes in TX amplitude and,
% when enabled, iteratively adjust the digital TX amplitude
% toward the target received signal peak.
% =========================================================================
if cfg.autocal_enable
    fprintf('\n--- RX LEVEL AUTOCALIBRATION ---\n');

    amp_original = cfg.tx_amplitude;


    [pk_a, ~] = medir_nivel(b, cfg, cfg.tx_amplitude * tone, 3);
    [pk_b, ~] = medir_nivel(b, cfg, (cfg.tx_amplitude/4) * tone, 3);
    cambio = abs(pk_a - pk_b) / max(pk_a, eps);
    fprintf('Response test: peak %.4f -> %.4f by dividing the amplitude by 4 (change %.1f%%)\n', ...
        pk_a, pk_b, 100*cambio);

    if cambio < 0.15
        warning(['The RX level does NOT correspond to the TX amplitude (change %.1f%%). ' ...
            'The receiver is not controlled by your transmitter. The ' ...
            ' self-calibration is aborted, and tx_amplitude is restored = %.4f.'], ...
            100*cambio, amp_original);
        cfg.tx_amplitude = amp_original;
        cfg.autocal_aborted = true;
    else
        cfg.autocal_aborted = false;
        for it = 1:cfg.autocal_max_iter
            [pk, ~] = medir_nivel(b, cfg, cfg.tx_amplitude * tone, cfg.autocal_blocks);
            fprintf('  iter %d: tx_amplitude=%.4f -> peak RX=%.3f\n', ...
                it, cfg.tx_amplitude, pk);

            if pk <= 0
                warning('RX peak is zero. Check wiring/antennas.');
                break;
            end
            if abs(pk - cfg.autocal_target_peak) < 0.08 * cfg.autocal_target_peak
                break;
            end

            factor = cfg.autocal_target_peak / pk;
            factor = max(0.25, min(4, factor));
            new_amp = cfg.tx_amplitude * factor;

            if new_amp > 0.5
                new_amp = 0.5;
                warning(['Digital amplitude at maximum (0.5) and the RX peak remains ' ...
                    'in %.3f. Increase cfg.rx_gain.'], pk);
            end
            if new_amp < 0.005
                new_amp = 0.005;
                warning(['Digital amplitude at its minimum usable level, and the RX peak remains ' ...
                    'in %.3f. Decrease cfg.rx_gain.'], pk);
            end

            cfg.tx_amplitude = new_amp;
            if new_amp >= 0.5 || new_amp <= 0.005
                break;
            end
        end
    end
    fprintf('  final amplitude: %.4f\n', cfg.tx_amplitude);
end

tx_block = cfg.tx_amplitude * tone;
fprintf('Warm-up complete. Starting real capture...\n');

% =========================================================================
% 5. CONTINUOUS I/Q ACQUISITION
% Transmit and receive data in consecutive blocks, append
% the available received samples, and record block-level
% measurements and streaming overruns.
% =========================================================================
total_overruns = 0;
fprintf('Audible indicator: The capture process begins.\n');
emitir_indicador_sonoro('inicio');
fprintf('Continuously capturing %d samples (%.1f s)...\n', ...
    total_samples, cfg.capture_seconds);

idx = 1;
try
    for k = 1:cfg.num_blocks
        b.transmit(tx_block, cfg.timeout_ms, 0, false, false);
        [rx, ~, actual_count, overrun] = b.receive(cfg.block_len);

        actual_counts(k) = actual_count;
        total_overruns = total_overruns + overrun;

        if actual_count > 0
            block_peak(k) = max(abs(rx(1:actual_count)));
            block_rms(k)  = sqrt(mean(abs(rx(1:actual_count)).^2));
        end

        if actual_count ~= cfg.block_len
            warning('Block %d is incomplete (%d/%d samples)', ...
                k, actual_count, cfg.block_len);
        end

        rx_continuous(idx:idx+actual_count-1) = rx(1:actual_count);
        idx = idx + actual_count;
    end
catch ME
    fprintf('\nERROR during capture\n');
    fprintf('%s\n', ME.message);
    rethrow(ME);
end

% =========================================================================
% 6. ACQUISITION INTEGRITY AND SIGNAL-LEVEL ASSESSMENT
% Evaluate the number of received samples, incomplete blocks,
% streaming overruns, and near-full-scale I/Q activity.
% =========================================================================
total_actual = idx - 1;
n_incomplete = nnz(actual_counts ~= cfg.block_len);
fprintf('\n--- TX/RX BLOCK REPORT ---\n');
fprintf('Incomplete blocks: %d of %d (%.2f%%).\n', ...
    n_incomplete, cfg.num_blocks, 100*n_incomplete/cfg.num_blocks);
fprintf('Actual samples captured: %d out of %d nominal (%.2f%%).\n', ...
    total_actual, total_samples, 100*total_actual/total_samples);
if total_actual < total_samples
    warning(['rx_continuous was truncated to %d samples (out of a nominal %d). ' ...
        'the time axis at a Fs=%.0f Hz no longer corresponds to the nominal %.1f s.'], ...
        total_actual, total_samples, cfg.fs, cfg.capture_seconds);
end
if total_actual == 0
    error('No valid sample was received.');
end
rx_continuous = rx_continuous(1:total_actual);

b.transmit(tx_zeros, cfg.timeout_ms, 0, false, true);

fprintf('\n--- CAPTURE REPORT ---\n');
if total_overruns > 0
    warning('DANGER: %d overruns. The phase has been disrupted. The radar will fail.', total_overruns);
else
    fprintf('Data continuity: EXCELLENT (0 overruns detected).\n');
end

idx_sat = abs(real(rx_continuous)) >= 0.99 | abs(imag(rx_continuous)) >= 0.99;
porcentaje = 100*nnz(idx_sat)/length(rx_continuous);
fprintf('Near-full-scale I/Q samples: %d (%.6f%%)\n', ...
    nnz(idx_sat), porcentaje);
if porcentaje > 0.01
    warning('Significant near-full-scale activity detected. Check RX gain and TX amplitude.');
end

pk_global = max(abs(rx_continuous));
fprintf('Global |rx| peak: %.3f (objective ~%.2f)\n', pk_global, cfg.autocal_target_peak);
if pk_global < 0.05
    warning(['Very low RX peak (%.3f): You are using less than 5%% of the ADC. Increase ' ...
        'cfg.rx_gain or cfg.tx_amplitude or the SNR will be terrible.'], pk_global);
end

fprintf('\n========== FULL BUFFER ==========\n');
fprintf('I: min %.4f   max %.4f\n', min(real(rx_continuous)), max(real(rx_continuous)));
fprintf('Q: min %.4f   max %.4f\n', min(imag(rx_continuous)), max(imag(rx_continuous)));

cfg.overruns = total_overruns;
cfg.saturation = porcentaje;
cfg.actual_counts = actual_counts;
cfg.block_peak = block_peak;
cfg.block_rms = block_rms;
cfg.n_incomplete_blocks = n_incomplete;
cfg.total_actual_samples = total_actual;
cfg.peak_global = pk_global;
cfg.capture_valid = (total_overruns == 0) && (total_actual > 0);

clear cleanup_obj;
emitir_indicador_sonoro('fin');
fprintf('Capture complete.\n');

% =========================================================================
% Auxiliar functions
% =========================================================================
% -------------------------------------------------------------------------
% MEDIR_NIVEL
% Measures the maximum complex-signal magnitude, maximum
% individual I/Q component, and RMS level over multiple
% received blocks.
% -------------------------------------------------------------------------
function [pk, rms_v, pk_comp] = medir_nivel(b, cfg, tx_blk, nblocks)
    pk = 0; pk_comp = 0; acc = 0; cnt = 0;
    for kk = 1:nblocks
        b.transmit(tx_blk, cfg.timeout_ms, 0, false, false);
        [rxm, ~, nc, ~] = b.receive(cfg.block_len);
        if nc > 0
            a = abs(rxm(1:nc));
            pk = max(pk, max(a));
            pk_comp = max(pk_comp, max([abs(real(rxm(1:nc))); abs(imag(rxm(1:nc)))]));
            acc = acc + sum(a.^2);
            cnt = cnt + nc;
        end
    end
    rms_v = sqrt(acc / max(cnt, 1));
end

% -------------------------------------------------------------------------
% VERIFICAR_GANANCIA_RX
% Compares the requested RX gain with the value reported
% by the SDR hardware and issues a warning or error
% according to the strict_gain configuration.
% -------------------------------------------------------------------------
function verificar_ganancia_rx(cfg)
    if ~isfinite(cfg.rx_gain_actual)
        return;
    end
    d = abs(cfg.rx_gain_actual - cfg.rx_gain);
    if d <= cfg.gain_tol_db
        return;
    end
    msg = sprintf([ ...
        'The requested RX gain was %.1f dB, but the hardware ' ...
        'reports %.1f dB after enabling RX.\n' ...
        'Check the AGC mode, the gain range, and the hardware response.'], ...
        cfg.rx_gain, cfg.rx_gain_actual);
    if cfg.strict_gain
        error('%s\n\n(To capture it as well, enter cfg.strict_gain = false.)', msg);
    else
        warning('%s', msg);
    end
end

% -------------------------------------------------------------------------
% CONFIGURE_BLADERF
% Configures the SDR RF channels, carrier frequencies,
% sampling rates, bandwidths, gain settings, AGC mode,
% and USB streaming buffers.
% -------------------------------------------------------------------------
function cfg = configure_bladerf(b, cfg)
    b.tx.channel = cfg.tx_channel;
    b.rx.channel = cfg.rx_channel;
    b.tx.frequency = cfg.fc;
    b.rx.frequency = cfg.fc;
    b.tx.samplerate = cfg.fs;
    b.rx.samplerate = cfg.fs;
    b.tx.bandwidth = cfg.bandwidth;
    b.rx.bandwidth = cfg.bandwidth;

    fprintf('\n--- GAIN MODE RX ---\n');

    b.rx.agc = cfg.rx_agc;

    modo_leido = b.rx.agc;

    fprintf('AGC RX requested: %s | read: %s\n', ...
        cfg.rx_agc, string(modo_leido));

    if ~strcmpi(string(modo_leido), string(cfg.rx_agc))
        error('Manual AGC mode could not be enabled.');
    end

    fprintf('\n--- GAINS ---\n');
    b.tx.gain = cfg.tx_gain;
    pause(0.05);
    try cfg.tx_gain_actual = b.tx.gain; catch, cfg.tx_gain_actual = NaN; end


    fprintf('TX requested %.2f dB -> real %.2f dB\n', ...
        cfg.tx_gain, cfg.tx_gain_actual);

    fprintf('RX requested %.2f dB -> would be applied when enabling RX\n', ...
        cfg.rx_gain);

    if isfinite(cfg.tx_gain_actual) && abs(cfg.tx_gain_actual - cfg.tx_gain) > 0.5
        warning(['The TX gain was %.2f dB instead of %.2f dB. Use ' ...
            'cfg.tx_amplitude to adjust the level in the digital domain.'], ...
            cfg.tx_gain_actual, cfg.tx_gain);
    end

    fprintf('-----------------\n');

    b.tx.config.num_buffers = cfg.num_buffers;
    b.tx.config.buffer_size = cfg.buffer_size;
    b.tx.config.num_transfers = cfg.num_transfers;
    b.tx.config.timeout_ms = cfg.timeout_ms;
    b.rx.config.num_buffers = cfg.num_buffers;
    b.rx.config.buffer_size = cfg.buffer_size;
    b.rx.config.num_transfers = cfg.num_transfers;
    b.rx.config.timeout_ms = cfg.timeout_ms;
end

% -------------------------------------------------------------------------
% EMITIR_INDICADOR_SONORO
% Generates an audible notification indicating the start
% or completion of the acquisition.
% -------------------------------------------------------------------------
function emitir_indicador_sonoro(tipo)
    try
        fs_audio = 8000;
        dur_tono = 0.12;
        dur_pausa = 0.06;
        t = (0:round(dur_tono*fs_audio)-1).' / fs_audio;
        envolvente = 10*sin(pi*(0:numel(t)-1).' / max(numel(t)-1, 1)).^2;
        silencio = zeros(round(dur_pausa*fs_audio), 1);

        if strcmpi(tipo, 'inicio')
            frecuencias = [700, 1000];
        else
            frecuencias = [1000, 700];
        end

        y = zeros(0,1);
        for f = frecuencias
            tono = envolvente .* sin(2*pi*f*t);
            y = [y; tono; silencio];
        end

        sound(y, fs_audio);
        pause(numel(y)/fs_audio + 0.15);
    catch ME
        warning('The audio indicator could not be played.: %s', ME.message);
    end
end

% -------------------------------------------------------------------------
% CLEANUP_BLADERF
% Stops TX/RX streaming and releases the SDR resources
% when the acquisition completes or an error occurs.
% -------------------------------------------------------------------------
function cleanup_bladerf(b, cfg)
    try
        b.transmit(complex(zeros(cfg.block_len, 1)), cfg.timeout_ms, 0, false, true);
    catch
    end
    try b.tx.stop(); catch, end
    pause(0.05);
    try b.rx.stop(); catch, end
    clear b;
end
end
