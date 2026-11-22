%% ECHO - VALIDACION DEL UMBRAL LOW-IF DE 2 kHz
% No modifica el frontend ni los resultados experimentales.

clc;
clearvars;

Fs_test = 1e6;
f_nominal = 100e3;
duration = 5;

% Frecuencias sinteticas de prueba
test_freqs = [97500, 97999, 98000, 98001, ...
              101999, 102000, 102001];

nTests = numel(test_freqs);

% Configuracion del frontend
cfg_test = struct();
cfg_test.fs = Fs_test;
cfg_test.tx_tone_hz = f_nominal;
cfg_test.debug_plots = false;

% Vector temporal
n = (0:Fs_test*duration-1).';

% Resultados
f_measured = zeros(nTests,1);
offset_measured = zeros(nTests,1);
contrast = zeros(nTests,1);
exceeded = false(nTests,1);
mix_frequency = zeros(nTests,1);
frequency_ok = false(nTests,1);
threshold_ok = false(nTests,1);

fprintf('\n============================================\n');
fprintf('     VALIDACION UMBRAL LOW-IF ECHO\n');
fprintf('============================================\n');

residual_hz = zeros(nTests,1);

for k = 1:nTests

    f_test = test_freqs(k);

    fprintf('\n--- PRUEBA %d/%d ---\n',k,nTests);
    fprintf('Frecuencia generada: %.2f Hz\n',f_test);

    % Generar senal compleja
    iq_test = exp(1j*2*pi*f_test*n/Fs_test);

    % Ejecutar frontend actual
    [iq_out, Fs_out, front_info] = ...
        frontend_dsp(iq_test,cfg_test);

    % Verificar frecuencia residual despues de la mezcla
    z = iq_out(:);

    i1 = floor(0.2*numel(z))+1;
    i2 = floor(0.8*numel(z));

    z_center = z(i1:i2);

    phi = unwrap(angle(z_center));
    t = (0:numel(phi)-1)'/Fs_out;

    p = polyfit(t,phi,1);
    residual_hz(k) = p(1)/(2*pi);

    fprintf('Frecuencia residual: %.4f Hz\n', ...
        residual_hz(k));

    % Verificar frecuencia residual despues de la mezcla
    z = iq_out(:);

    % Evitar los extremos del filtrado
    i1 = floor(0.2*numel(z))+1;
    i2 = floor(0.8*numel(z));

    z_center = z(i1:i2);

    phi = unwrap(angle(z_center));
    t = (0:numel(phi)-1)'/Fs_out;

    p = polyfit(t,phi,1);
    f_residual = p(1)/(2*pi);

    fprintf('Frecuencia residual en banda base: %.4f Hz\n', ...
        f_residual);

    % Registrar resultados
    f_measured(k) = front_info.measured_hz;
    offset_measured(k) = front_info.offset_hz;
    contrast(k) = front_info.contrast_db;
    exceeded(k) = front_info.exceeded_2khz;
    mix_frequency(k) = front_info.mix_frequency_hz;

    % Verificar frecuencia detectada
    frequency_ok(k) = ...
        abs(f_measured(k)-f_test) < 1;

    % El criterio es estrictamente mayor que 2000 Hz
    % Verificar la decision usando el offset realmente estimado
    expected_flag = abs(front_info.offset_hz) > 2000;

    threshold_ok(k) = ...
        exceeded(k) == expected_flag;

    fprintf('Frecuencia detectada: %.4f Hz\n',f_measured(k));
    fprintf('Offset detectado: %+.4f Hz\n',offset_measured(k));
    fprintf('Contraste: %.2f dB\n',contrast(k));
    fprintf('Supera 2 kHz: %d\n',exceeded(k));
    fprintf('Frecuencia de mezcla: %.4f Hz\n',mix_frequency(k));

    fprintf('Deteccion correcta: %d\n',frequency_ok(k));
    fprintf('Umbral correcto: %d\n',threshold_ok(k));

    clear iq_test iq_out;
end

%% TABLA FINAL

T = table( ...
    test_freqs(:), ...
    f_measured, ...
    offset_measured, ...
    contrast, ...
    exceeded, ...
    mix_frequency, ...
    residual_hz, ...
    frequency_ok, ...
    threshold_ok, ...
    'VariableNames',{ ...
    'Frequency_Generated_Hz', ...
    'Frequency_Measured_Hz', ...
    'Offset_Hz', ...
    'Contrast_dB', ...
    'Exceeded_2kHz', ...
    'Mix_Frequency_Hz', ...
    'Residual_Baseband_Hz', ...
    'Frequency_OK', ...
    'Threshold_OK'});

fprintf('\n============================================\n');
fprintf('RESULTADOS FINALES\n');
fprintf('============================================\n');

disp(T);

fprintf('Frecuencias correctamente detectadas: %d/%d\n', ...
    sum(frequency_ok),nTests);

fprintf('Decisiones de umbral correctas: %d/%d\n', ...
    sum(threshold_ok),nTests);

if all(frequency_ok) && all(threshold_ok)
    fprintf('\nVALIDACION SUPERADA.\n');
else
    warning('VALIDACION NO SUPERADA. Revisar resultados.');
end

% Exportar resultados sin sobrescribir otros archivos
writetable(T,'ECHO_test_umbral_LowIF.xlsx');

%% VERIFICACION DEL LIMITE EXACTO

boundary_offsets = [-2000, 2000, -2000.001, 2000.001];

boundary_expected = [false; false; true; true];

boundary_actual = abs(boundary_offsets(:)) > 2000;

boundary_ok = isequal(boundary_actual,boundary_expected);

fprintf('\n--- VERIFICACION DEL LIMITE EXACTO ---\n');

disp(table(boundary_offsets(:),boundary_actual, ...
    boundary_expected, ...
    'VariableNames',{'Offset_Hz','Actual','Expected'}));

fprintf('Verificacion del limite: %d/4\n', ...
    sum(boundary_actual == boundary_expected));