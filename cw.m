function [rx_continuous,cfg] = cw()
% =========================================================================
% CAPTURA ROBUSTA BLADERF - RADAR CW (ONDA CONTINUA)
% =========================================================================

cfg = struct();
cfg.serial = '8c14dec0ebc947b88424691d0dbf744f'; 

cfg.fc = 5.8e9;              % Frecuencia portadora (Hz)
cfg.fs = 1e6;                % Tasa de muestreo (Sa/s)
cfg.bandwidth = 1e6;         % Ancho de banda (Hz)
cfg.tx_tone_hz = 100e3;      % Tono Low-IF para evitar ruido 1/f
cfg.tx_amplitude = 0.15;     % Amplitud (0 - 0.5)

cfg.tx_channel = 'TX1';
cfg.rx_channel = 'RX2';      
cfg.tx_gain = -45;           % dB
cfg.rx_gain = -15;           % dB (bajado de -5: el pico por bloque estaba
                              % pegado a 1.0 en toda la captura de 40s)
cfg.rx_agc = 'manual';

cfg.block_len = 10000;
cfg.capture_seconds = 40;
cfg.num_blocks = ceil(cfg.capture_seconds * cfg.fs / cfg.block_len);
cfg.capture_seconds = (cfg.block_len * cfg.num_blocks) / cfg.fs;
cfg.num_buffers = 32;
cfg.buffer_size = 30720;
cfg.num_transfers = 8;
cfg.timeout_ms = 5000;

assert(cfg.tx_amplitude > 0 && cfg.tx_amplitude <= 0.5, ...
    'cfg.tx_amplitude debe estar entre 0 y 0.5.');

samples_per_period = cfg.fs / cfg.tx_tone_hz;
assert(mod(cfg.block_len, samples_per_period) == 0, ...
    'cfg.block_len debe ser múltiplo entero de fs/tx_tone_hz para fase continua entre bloques.');


if exist('bladeRF', 'class') ~= 8 && exist('bladeRF', 'file') ~= 2
    error('No se encontro bladeRF.m en el path.');
end

fprintf('Abriendo bladeRF...\n');
b = bladeRF(['*:serial=' cfg.serial]);
cleanup_obj = onCleanup(@() cleanup_bladerf(b, cfg));

configure_bladerf(b, cfg);

total_samples = cfg.block_len * cfg.num_blocks;
rx_continuous = complex(zeros(total_samples, 1)); 
actual_counts = zeros(cfg.num_blocks, 1);
% Envolvente de amplitud por bloque, para diagnosticar la captura completa
% (no solo los primeros ms). Se llenan en el loop principal, un valor por
% cada intento de transmit+receive, independientemente de si el bloque
% resulto completo o no.
block_peak = zeros(cfg.num_blocks, 1);
block_rms  = zeros(cfg.num_blocks, 1);
% --- CÓDIGO CORREGIDO ---
%% Generar bloque TX primero, en la tranquilidad de la RAM
n = (0:cfg.block_len-1).';
tx_block = cfg.tx_amplitude * exp(1j*2*pi*cfg.tx_tone_hz*n/cfg.fs);

% =========================================================================
% 1. WARM-UP MIENTRAS EL PUERTO ESTÁ CERRADO Y SEGURO
% =========================================================================
fprintf('Esperando estabilización de osciladores a 5.8 GHz (Cold Start)...\n');
% La pausa ocurre aquí, el hardware se calienta, pero el USB no envía nada aún
pause(1.0); 

% Ahora sí, encender y alimentar el hardware inmediatamente sin pausas
fprintf('Iniciando RX %s y TX %s...\n', cfg.rx_channel, cfg.tx_channel);
b.rx.start();
b.tx.start();

fprintf('Limpiando buffers USB...\n');
% 2. Descartar bloques iniciales (Warm-up de buffers)
for k = 1:30
    b.transmit(tx_block, cfg.timeout_ms, 0, false, false);
    [~, ~, ~, ov_warmup] = b.receive(cfg.block_len);
    if ov_warmup
        warning('Overrun durante warm-up en iteración %d', k);
    end
end
fprintf('Calentamiento listo. Iniciando captura real...\n');


total_overruns = 0;
% actual_counts = zeros(cfg.num_blocks, 1);

fprintf('Capturando %d muestras de forma continua...\n', total_samples);

idx = 1;

try
    for k = 1:cfg.num_blocks
        b.transmit(tx_block, cfg.timeout_ms, 0, false, false);
        [rx, ~, actual_count, overrun] = b.receive(cfg.block_len);


        actual_counts(k) = actual_count;
        total_overruns = total_overruns + overrun;

        % Envolvente real de este bloque, tal como llego del hardware.
        % Se calcula ANTES de cualquier recorte, para que un bloque corto
        % o vacio quede reflejado como pico/RMS bajo (no como NaN u
        % omision silenciosa).
        if actual_count > 0
            block_peak(k) = max(abs(rx(1:actual_count)));
            block_rms(k)  = sqrt(mean(abs(rx(1:actual_count)).^2));
        end

        % Verificar cantidad recibida
        if actual_count ~= cfg.block_len
            warning("Bloque %d incompleto (%d/%d muestras)",...
                k,actual_count,cfg.block_len);
        end

        rx_continuous(idx:idx+actual_count-1)=rx(1:actual_count);

        idx=idx+actual_count;

    end
catch ME

    fprintf('\nERROR durante la captura\n');
    fprintf('%s\n',ME.message);

    rethrow(ME);

end

% Recortar el buffer al total de muestras REALMENTE recibidas. Si algun
% bloque salio incompleto, el resto del buffer preasignado se habria
% quedado en ceros literales al final del vector -- eso diluye el
% reporte de saturacion (divide entre una longitud inflada) y descuadra
% el eje de tiempo asumido (Fs constante) frente al tiempo real
% transcurrido, sin avisar en ningun lado.
total_actual = idx - 1;
n_incomplete = nnz(actual_counts ~= cfg.block_len);
fprintf('\n--- REPORTE DE BLOQUES TX/RX ---\n');
fprintf('Bloques incompletos: %d de %d (%.2f%%).\n', ...
    n_incomplete, cfg.num_blocks, 100*n_incomplete/cfg.num_blocks);
fprintf('Muestras reales capturadas: %d de %d nominales (%.2f%%).\n', ...
    total_actual, total_samples, 100*total_actual/total_samples);
if total_actual < total_samples
    warning(['rx_continuous se recorto a %d muestras (de %d nominales) por ' ...
        'bloques incompletos. El eje de tiempo asumido a Fs=%.0f Hz ya NO ' ...
        'corresponde exactamente a los %.1f s nominales de principio a fin.'], ...
        total_actual, total_samples, cfg.fs, cfg.capture_seconds);
end
rx_continuous = rx_continuous(1:max(total_actual, 1));

% Detener transmisión de forma segura
b.transmit(complex(zeros(cfg.block_len, 1)), cfg.timeout_ms, 0, false, true);


% =========================================================================
% 3. VERIFICACIÓN DE ROBUSTEZ Y CALIDAD DE SEÑAL
% =========================================================================
fprintf('\n--- REPORTE DE CAPTURA ---\n');

if total_overruns > 0
    warning('PELIGRO: Se detectaron %d overruns. La fase se ha roto. El radar fallará.', total_overruns);
else
    fprintf('Continuidad de datos: EXCELENTE (0 overruns detectados).\n');
end


% =========================================================================
% REPORTE DE SATURACIÓN
% =========================================================================

idx_sat = abs(real(rx_continuous))>=0.99 | ...
          abs(imag(rx_continuous))>=0.99;

porcentaje = 100*nnz(idx_sat)/length(rx_continuous);

fprintf("Muestras saturadas: %d (%.6f%%)\n",...
    nnz(idx_sat),porcentaje);

if porcentaje > 0.01
    warning("Existe saturación importante.");
end

fprintf("\n========== BUFFER COMPLETO ==========\n");

fprintf("I: min %.4f   max %.4f\n", ...
    min(real(rx_continuous)), ...
    max(real(rx_continuous)));

fprintf("Q: min %.4f   max %.4f\n", ...
    min(imag(rx_continuous)), ...
    max(imag(rx_continuous)));
% AHORA SÍ: Removemos el DC residual (Fuga de hardware) ANTES de la demodulación
rx_continuous = rx_continuous - mean(rx_continuous);

cfg.overruns = total_overruns;
cfg.saturation = porcentaje;
cfg.actual_counts = actual_counts;
cfg.block_peak = block_peak;
cfg.block_rms = block_rms;
cfg.n_incomplete_blocks = n_incomplete;
cfg.total_actual_samples = total_actual;

% Limpieza final
clear cleanup_obj;
fprintf('Captura finalizada.\n');

% =========================================================================
% FUNCIONES AUXILIARES
% =========================================================================
function configure_bladerf(b, cfg)
    b.tx.channel = cfg.tx_channel;
    b.rx.channel = cfg.rx_channel;
    b.tx.frequency = cfg.fc;
    b.rx.frequency = cfg.fc;
    b.tx.samplerate = cfg.fs;
    b.rx.samplerate = cfg.fs;
    b.tx.bandwidth = cfg.bandwidth;
    b.rx.bandwidth = cfg.bandwidth;
    b.tx.gain = cfg.tx_gain;
    b.rx.agc = cfg.rx_agc;
    b.rx.gain = cfg.rx_gain;
    b.tx.config.num_buffers = cfg.num_buffers;
    b.tx.config.buffer_size = cfg.buffer_size;
    b.tx.config.num_transfers = cfg.num_transfers;
    b.tx.config.timeout_ms = cfg.timeout_ms;
    b.rx.config.num_buffers = cfg.num_buffers;
    b.rx.config.buffer_size = cfg.buffer_size;
    b.rx.config.num_transfers = cfg.num_transfers;
    b.rx.config.timeout_ms = cfg.timeout_ms;
    b.tx.channel = cfg.tx_channel;
    b.rx.channel = cfg.rx_channel;
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
