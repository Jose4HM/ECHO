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
cfg.tx_amplitude = 0.07;     % Amplitud (0 - 0.5)

cfg.tx_channel = 'TX1';
cfg.rx_channel = 'RX2';      
cfg.tx_gain = -45;           % dB
cfg.rx_gain = -15;            % dB 
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
    b.transmit(tx_block);
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
        % transmit() sin flags SOB/EOB evita cerrar el burst en cada
        % bloque, reduciendo discontinuidades periodicas en fase.
        b.transmit(tx_block);
        [rx, ~, actual_count, overrun] = b.receive(cfg.block_len);


        actual_counts(k) = actual_count;
        total_overruns = total_overruns + overrun;

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

% Detener transmisión de forma segura
if idx <= total_samples
    rx_continuous = rx_continuous(1:idx-1);
end

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
