function [rx_continuous, cfg] = cw()
%CW Captura continua para radar CW con bladeRF 2.0.
%   Genera un tono Low-IF, captura I/Q crudo y registra diagnosticos de
%   continuidad, saturacion y nivel antes del procesamiento fisiologico.

cfg = struct();
cfg.serial = '8c14dec0ebc947b88424691d0dbf744f';

cfg.fc = 5.94e9;             % Portadora (Hz). 5.8 GHz estaba inutilizable:
                             % el barrido de banda mostro WiFi saturando en
                             % 5.79 / 5.82 / 5.85 GHz (scan_banda.m)
cfg.fs = 1e6;                % Tasa de muestreo (Sa/s)
cfg.bandwidth = 1e6;         % Ancho de banda (Hz)
cfg.tx_tone_hz = 100e3;      % Tono Low-IF para evitar ruido 1/f
cfg.tx_amplitude = 0.25;     % Amplitud digital (0 - 0.5)

cfg.tx_channel = 'TX1';
cfg.rx_channel = 'RX2';
cfg.tx_gain = 65;            % dB. El TX es el UNICO mando de nivel util en
                             % esta placa. A -24 dB el eco del pecho
                             % quedaba enterrado en el ruido, solo la
                             % respiración podía ser hallada
cfg.rx_gain = 0;             % dB. Sin efecto en este binding, no aplica.
cfg.rx_agc  = 'manual';

% la escritura de ganancia RX no llega al
% AD9361 aunque se lea de vuelta bien. Se comprobo tres veces, incluida
% una con WiFi entrando a 0.4 de fondo de escala: el nivel no cambia al
% mover rx_gain 72 dB. Por tanto el guardian estricto solo estorba aqui.
% Poner true solo si el binding aplica la ganancia RX de forma fiable.
cfg.strict_gain = false;
cfg.gain_tol_db = 3;

% Linea base con TX apagado
cfg.baseline_blocks = 10;
% Tiempo con el tono TX activo antes de guardar datos. Permite que la fase
% RF/LO y la cola de streaming se estabilicen fuera de la captura util.
cfg.tx_settle_seconds = 20;

cfg.autocal_enable = false;
cfg.autocal_target_peak = 0.45;
cfg.autocal_max_iter = 6;
cfg.autocal_blocks = 5;

cfg.block_len = 10000;
cfg.capture_seconds = 60;
cfg.num_blocks = ceil(cfg.capture_seconds * cfg.fs / cfg.block_len);
cfg.capture_seconds = (cfg.block_len * cfg.num_blocks) / cfg.fs;
cfg.num_buffers = 32;
cfg.buffer_size = 30720;
cfg.num_transfers = 8;
cfg.timeout_ms = 5000;

assert(cfg.tx_amplitude > 0 && cfg.tx_amplitude <= 0.5, ...
    'cfg.tx_amplitude debe estar entre 0 y 0.5.');

samples_per_period = cfg.fs / cfg.tx_tone_hz;
assert(abs(samples_per_period - round(samples_per_period)) < 1e-12, ...
    'fs/tx_tone_hz debe ser entero para repetir el tono sin salto de fase.');
samples_per_period = round(samples_per_period);
assert(mod(cfg.block_len, samples_per_period) == 0, ...
    'cfg.block_len debe contener un numero entero de periodos del tono.');

if exist('bladeRF', 'class') ~= 8 && exist('bladeRF', 'file') ~= 2
    error('No se encontro bladeRF.m en el path.');
end

fprintf('Abriendo bladeRF...\n');
b = bladeRF(['*:serial=' cfg.serial]);
cleanup_obj = onCleanup(@() cleanup_bladerf(b, cfg));

cfg = configure_bladerf(b, cfg);

total_samples = cfg.block_len * cfg.num_blocks;
rx_continuous = complex(zeros(total_samples, 1));
actual_counts = zeros(cfg.num_blocks, 1);
block_peak = zeros(cfg.num_blocks, 1);
block_rms  = zeros(cfg.num_blocks, 1);

n = (0:cfg.block_len-1).';
tone = exp(1j*2*pi*cfg.tx_tone_hz*n/cfg.fs);
tx_block = cfg.tx_amplitude * tone;
tx_zeros = complex(zeros(cfg.block_len, 1));

% =========================================================================
% 1. WARM-UP
% =========================================================================
fprintf('Esperando bloqueo inicial de sintetizadores a %.3f GHz...\n', cfg.fc/1e9);
pause(1.0);

fprintf('Iniciando RX %s y TX %s...\n', cfg.rx_channel, cfg.tx_channel);
b.rx.start();
b.tx.start();

% Reafirmar la ganancia despues de arrancar el stream.
% cfg = reafirmar_ganancia_rx(b, cfg);

fprintf('Limpiando buffers USB...\n');
for k = 1:30
    b.transmit(tx_zeros, cfg.timeout_ms, 0, false, false);
    [~, ~, ~, ov_warmup] = b.receive(cfg.block_len);
    if ov_warmup
        warning('Overrun durante warm-up en iteracion %d', k);
    end
end

% =========================================================================
% 2. LINEA BASE CON TX APAGADO
% =========================================================================
fprintf('\n--- LINEA BASE (TX APAGADO) ---\n');
[pk_off, rms_off, pk_comp_off] = medir_nivel(b, cfg, tx_zeros, cfg.baseline_blocks);
fprintf('Pico |IQ| sin TX: %.4f   pico componente I/Q: %.4f\n', pk_off, pk_comp_off);
fprintf('RMS   sin TX: %.5f (%.1f dBFS)\n', rms_off, 20*log10(rms_off + 1e-12));

% Al pasar de ceros al tono, la cola TX y la fase RF necesitan tiempo para
% estabilizarse. Todo lo recibido aqui se descarta deliberadamente.
fprintf('\n--- ESTABILIZACION TX ENCENDIDO (%.1f s) ---\n', cfg.tx_settle_seconds);
n_settle = ceil(cfg.tx_settle_seconds * cfg.fs / cfg.block_len);
settle_overruns = 0;
settle_incomplete = 0;
for k = 1:n_settle
    b.transmit(tx_block, cfg.timeout_ms, 0, false, false);
    [~, ~, nc_settle, ov_settle] = b.receive(cfg.block_len);
    settle_overruns = settle_overruns + double(ov_settle ~= 0);
    settle_incomplete = settle_incomplete + double(nc_settle ~= cfg.block_len);
end
cfg.tx_settle_actual_seconds = n_settle * cfg.block_len / cfg.fs;
cfg.tx_settle_overruns = settle_overruns;
cfg.tx_settle_incomplete_blocks = settle_incomplete;
fprintf('Estabilizacion terminada: %.1f s, overruns %d, bloques incompletos %d.\n', ...
    cfg.tx_settle_actual_seconds, settle_overruns, settle_incomplete);
if settle_overruns > 0 || settle_incomplete > 0
    warning(['Hubo incidencias durante la estabilizacion TX. No pertenecen a ' ...
        'la captura util, pero conviene revisar el streaming.']);
end

[pk_on, rms_on, pk_comp_on] = medir_nivel(b, cfg, tx_block, cfg.baseline_blocks);
fprintf('Pico |IQ| con TX: %.4f   pico componente I/Q: %.4f\n', pk_on, pk_comp_on);
fprintf('RMS   con TX: %.5f (%.1f dBFS)\n', rms_on, 20*log10(rms_on + 1e-12));
razon_db = 20*log10((rms_on + 1e-12) / (rms_off + 1e-12));
fprintf('Incremento RMS al activar TX: %.1f dB\n', razon_db);

cfg.baseline_peak_off = pk_off;
cfg.baseline_rms_off  = rms_off;
cfg.baseline_peak_on  = pk_on;
cfg.baseline_peak_component_off = pk_comp_off;
cfg.baseline_peak_component_on  = pk_comp_on;
cfg.baseline_rms_on   = rms_on;
cfg.leak_over_floor_db = razon_db;

if pk_comp_off >= 0.99
    error(['El receptor satura con TX digitalmente apagado (pico I/Q %.3f). ' ...
        'Revise ganancia RX e interferencia externa cerca de %.3f GHz antes de capturar.'], ...
        pk_comp_off, cfg.fc/1e9);
end

if razon_db < 6
    warning(['El RMS solo aumenta %.1f dB al activar TX. Verifique que el ' ...
        'tono Low-IF este presente y que TX/antenas funcionen correctamente.'], razon_db);
end

% =========================================================================
% 3. AUTOCALIBRACION DE NIVEL (con verificacion de respuesta)
% =========================================================================
if cfg.autocal_enable
    fprintf('\n--- AUTOCALIBRACION DE NIVEL RX ---\n');

    amp_original = cfg.tx_amplitude;


    [pk_a, ~] = medir_nivel(b, cfg, cfg.tx_amplitude * tone, 3);
    [pk_b, ~] = medir_nivel(b, cfg, (cfg.tx_amplitude/4) * tone, 3);
    cambio = abs(pk_a - pk_b) / max(pk_a, eps);
    fprintf('Prueba de respuesta: pico %.4f -> %.4f al dividir amplitud /4 (cambio %.1f%%)\n', ...
        pk_a, pk_b, 100*cambio);

    if cambio < 0.15
        warning(['El nivel RX NO responde a la amplitud de TX (cambio %.1f%%). ' ...
            'Lo que domina el receptor no es tu transmisor. Se aborta la ' ...
            'autocalibracion y se restaura tx_amplitude = %.4f.'], ...
            100*cambio, amp_original);
        cfg.tx_amplitude = amp_original;
        cfg.autocal_aborted = true;
    else
        cfg.autocal_aborted = false;
        for it = 1:cfg.autocal_max_iter
            [pk, ~] = medir_nivel(b, cfg, cfg.tx_amplitude * tone, cfg.autocal_blocks);
            fprintf('  iter %d: tx_amplitude=%.4f -> pico RX=%.3f\n', ...
                it, cfg.tx_amplitude, pk);

            if pk <= 0
                warning('Pico RX nulo. Revisa cableado/antenas.');
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
                warning(['Amplitud digital al maximo (0.5) y el pico RX sigue ' ...
                    'en %.3f. Sube cfg.rx_gain.'], pk);
            end
            if new_amp < 0.005
                new_amp = 0.005;
                warning(['Amplitud digital en su minimo util y el pico RX sigue ' ...
                    'en %.3f. Baja cfg.rx_gain.'], pk);
            end

            cfg.tx_amplitude = new_amp;
            if new_amp >= 0.5 || new_amp <= 0.005
                break;
            end
        end
    end
    fprintf('  amplitud final: %.4f\n', cfg.tx_amplitude);
end

tx_block = cfg.tx_amplitude * tone;
fprintf('Calentamiento listo. Iniciando captura real...\n');

% =========================================================================
% 4. CAPTURA
% =========================================================================
total_overruns = 0;
fprintf('Indicador sonoro: comienza la captura util.\n');
emitir_indicador_sonoro('inicio');
fprintf('Capturando %d muestras de forma continua (%.1f s)...\n', ...
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
            warning('Bloque %d incompleto (%d/%d muestras)', ...
                k, actual_count, cfg.block_len);
        end

        rx_continuous(idx:idx+actual_count-1) = rx(1:actual_count);
        idx = idx + actual_count;
    end
catch ME
    fprintf('\nERROR durante la captura\n');
    fprintf('%s\n', ME.message);
    rethrow(ME);
end

total_actual = idx - 1;
n_incomplete = nnz(actual_counts ~= cfg.block_len);
fprintf('\n--- REPORTE DE BLOQUES TX/RX ---\n');
fprintf('Bloques incompletos: %d de %d (%.2f%%).\n', ...
    n_incomplete, cfg.num_blocks, 100*n_incomplete/cfg.num_blocks);
fprintf('Muestras reales capturadas: %d de %d nominales (%.2f%%).\n', ...
    total_actual, total_samples, 100*total_actual/total_samples);
if total_actual < total_samples
    warning(['rx_continuous se recorto a %d muestras (de %d nominales). El ' ...
        'eje de tiempo a Fs=%.0f Hz ya no corresponde a los %.1f s nominales.'], ...
        total_actual, total_samples, cfg.fs, cfg.capture_seconds);
end
if total_actual == 0
    error('No se recibio ninguna muestra valida.');
end
rx_continuous = rx_continuous(1:total_actual);

b.transmit(tx_zeros, cfg.timeout_ms, 0, false, true);

% =========================================================================
% 5. REPORTE DE CALIDAD
% =========================================================================
fprintf('\n--- REPORTE DE CAPTURA ---\n');
if total_overruns > 0
    warning('PELIGRO: %d overruns. La fase se rompio. El radar fallara.', total_overruns);
else
    fprintf('Continuidad de datos: EXCELENTE (0 overruns detectados).\n');
end

idx_sat = abs(real(rx_continuous)) >= 0.99 | abs(imag(rx_continuous)) >= 0.99;
porcentaje = 100*nnz(idx_sat)/length(rx_continuous);
fprintf('Muestras saturadas: %d (%.6f%%)\n', nnz(idx_sat), porcentaje);
if porcentaje > 0.01
    warning('Saturacion importante. Baja cfg.rx_gain o cfg.tx_amplitude.');
end

pk_global = max(abs(rx_continuous));
fprintf('Pico |rx| global: %.3f (objetivo ~%.2f)\n', pk_global, cfg.autocal_target_peak);
if pk_global < 0.05
    warning(['Pico RX muy bajo (%.3f): usas menos del 5%% del ADC. Sube ' ...
        'cfg.rx_gain o cfg.tx_amplitude o el SNR sera pesimo.'], pk_global);
end

fprintf('\n========== BUFFER COMPLETO ==========\n');
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
fprintf('Captura finalizada.\n');

% =========================================================================
% FUNCIONES AUXILIARES
% =========================================================================
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

function cfg = reafirmar_ganancia_rx(b, cfg)
    try
        % b.rx.gain = cfg.rx_gain;
        pause(0.1);
        g = b.rx.gain;
    catch
        return;
    end
    if abs(g - cfg.rx_gain_actual) > 0.5
        fprintf('RX gain tras arrancar el stream: %.2f dB (antes %.2f dB)\n', ...
            g, cfg.rx_gain_actual);
    end
    cfg.rx_gain_actual = g;
    verificar_ganancia_rx(cfg);
end

function verificar_ganancia_rx(cfg)
    if ~isfinite(cfg.rx_gain_actual)
        return;
    end
    d = abs(cfg.rx_gain_actual - cfg.rx_gain);
    if d <= cfg.gain_tol_db
        return;
    end
    msg = sprintf([ ...
        'La ganancia RX NO se aplico: pediste %.1f dB y el hardware reporta ' ...
        '%.1f dB.\nEso no es un recorte a rango, es que la escritura no tomo ' ...
        'efecto.\nCausa mas probable: el AGC sigue activo, o el binding usa ' ...
        'otro nombre de propiedad.\nCorre  diag_ganancia(''%s'')  para ver que ' ...
        'expone tu version.'], cfg.rx_gain, cfg.rx_gain_actual, cfg.serial);
    if cfg.strict_gain
        error('%s\n\n(Para capturar igualmente, pon cfg.strict_gain = false.)', msg);
    else
        warning('%s', msg);
    end
end

function cfg = configure_bladerf(b, cfg)
    b.tx.channel = cfg.tx_channel;
    b.rx.channel = cfg.rx_channel;
    b.tx.frequency = cfg.fc;
    b.rx.frequency = cfg.fc;
    b.tx.samplerate = cfg.fs;
    b.rx.samplerate = cfg.fs;
    b.tx.bandwidth = cfg.bandwidth;
    b.rx.bandwidth = cfg.bandwidth;

    fprintf('\n--- MODO DE GANANCIA RX ---\n');
    modo_puesto = '';
    for nombre = {'agc', 'gain_mode'}
        campo = nombre{1};
        try
            b.rx.(campo) = cfg.rx_agc;
            leido = b.rx.(campo);
            fprintf('b.rx.%-10s = ''%s'' -> lee ''%s''\n', campo, cfg.rx_agc, string(leido));
            modo_puesto = campo;
        catch
            fprintf('b.rx.%-10s no disponible en este binding\n', campo);
        end
    end
    cfg.rx_gain_mode_prop = modo_puesto;
    if isempty(modo_puesto)
        warning(['No se pudo fijar el modo de ganancia manual. Si el AGC ' ...
            'queda activo, cfg.rx_gain no tendra ningun efecto.']);
    end

    fprintf('\n--- GANANCIAS ---\n');
    b.tx.gain = cfg.tx_gain;
    pause(0.05);
    try cfg.tx_gain_actual = b.tx.gain; catch, cfg.tx_gain_actual = NaN; end

    b.rx.gain = cfg.rx_gain;
    pause(0.1);
    try g = b.rx.gain; catch, g = NaN; end
    if isfinite(g) && abs(g - cfg.rx_gain) > cfg.gain_tol_db
        fprintf('Primer intento fallido (%.1f dB). Reintentando...\n', g);
        try
            if ~isempty(modo_puesto)
                b.rx.(modo_puesto) = cfg.rx_agc;
            end
            pause(0.1);
            b.rx.gain = cfg.rx_gain;
            pause(0.1);
            g = b.rx.gain;
        catch
        end
    end
    cfg.rx_gain_actual = g;

    fprintf('TX solicitado %.2f dB -> real %.2f dB\n', cfg.tx_gain, cfg.tx_gain_actual);
    fprintf('RX solicitado %.2f dB -> real %.2f dB\n', cfg.rx_gain, cfg.rx_gain_actual);
    if isfinite(cfg.tx_gain_actual) && abs(cfg.tx_gain_actual - cfg.tx_gain) > 0.5
        warning(['La ganancia TX quedo en %.2f dB en vez de %.2f dB. Usa ' ...
            'cfg.tx_amplitude para ajustar el nivel en el dominio digital.'], ...
            cfg.tx_gain_actual, cfg.tx_gain);
    end
    verificar_ganancia_rx(cfg);
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

function emitir_indicador_sonoro(tipo)
    % Dos tonos breves: ascendentes al iniciar y descendentes al terminar.
    % La pausa espera a que acabe el audio para no incluirlo en la captura.
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
            y = [y; tono; silencio]; %#ok<AGROW>
        end

        sound(y, fs_audio);
        pause(numel(y)/fs_audio + 0.15);
    catch ME
        warning('No se pudo reproducir el indicador sonoro: %s', ME.message);
    end
end

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
