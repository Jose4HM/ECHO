%% test_sintetico_radar.m
% =========================================================================
% BANCO DE PRUEBAS SINTETICO PARA RADAR CW DE SIGNOS VITALES
%
% No requiere bladeRF.
%
% Prueba dos niveles de la cadena:
%   A) Front-end Low-IF -> banda base (opcional, corto para no usar mucha RAM)
%   B) Procesamiento fisiologico completo a 100 Hz con varios casos conocidos
%
% Modelo de banda base:
%
% z(t) = exp(j*phi_LO(t)) * ...
%        [ C + A*exp(j*4*pi/lambda*(x_resp(t)+x_hr(t)+x_mov(t))) ] + ruido
%
% Esto reproduce:
%   - fuga TX->RX / clutter estatico (C)
%   - arco IQ del eco del pecho (A)
%   - respiracion
%   - latido
%   - offset y deriva entre sintetizadores
%   - ruido complejo
%   - movimiento corporal puntual
%
% Requisitos:
%   frontend_dsp.m
%   procesar_signos_vitales_fixed.m
%   Signal Processing Toolbox
%
% IMPORTANTE:
% Los nombres de archivo deben coincidir con los nombres de las funciones.
% No uses "frontend_dsp(1).m" ni "procesar_signos_vitales_fixed(2).m".
% Renombralos a:
%   frontend_dsp.m
%   procesar_signos_vitales_fixed.m
% =========================================================================

clear; clc; close all;
rng(42);                         % prueba reproducible

%% CONFIGURACION GENERAL
fc = 5.94e9;
c0 = 299792458;
lambda = c0 / fc;

Fs_bb = 100;                     % salida esperada del frontend
duracion_s = 30;

% Mismos filtros fisiologicos previstos para el sistema final.
fc_hp = 0.04;
dHP = designfilt('highpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', fc_hp, 'SampleRate', Fs_bb);

dLP_rr = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 0.8, 'SampleRate', Fs_bb);

fprintf('============================================================\n');
fprintf(' BANCO DE PRUEBAS SINTETICO - RADAR CW\n');
fprintf(' fc = %.3f GHz | lambda = %.2f mm | Fs BB = %.0f Hz\n', ...
    fc/1e9, lambda*1e3, Fs_bb);
fprintf('============================================================\n\n');

%% =========================================================================
% A) PRUEBA CORTA DEL FRONT-END
% =========================================================================
PROBAR_FRONTEND = true;

if PROBAR_FRONTEND
    if exist('frontend_dsp', 'file') ~= 2
        warning(['No se encontro frontend_dsp.m. Se omite la prueba del ' ...
            'front-end. Revisa que el archivo tenga exactamente ese nombre.']);
    else
        fprintf('--- A) PRUEBA DEL FRONT-END LOW-IF ---\n');

        cfg = struct();
        cfg.fs = 1e6;
        cfg.tx_tone_hz = 100e3;
        cfg.debug_plots = false;

        % Solo 3 s para limitar memoria. Esta prueba NO estima RPM/BPM:
        % valida mezclador, estimacion de tono y decimacion.
        dur_fe = 3;
        Nfe = round(dur_fe * cfg.fs);
        tfe = (0:Nfe-1).' / cfg.fs;

        tone_offset_hz = 23.0;
        f_mod = 1.20;            % modulacion conocida = 72 ciclos/min
        phase_mod = 0.35 * sin(2*pi*f_mod*tfe);

        Cfe = 0.25 * exp(1j*0.25);
        Afe = 0.040;

        raw = exp(1j*2*pi*(cfg.tx_tone_hz + tone_offset_hz)*tfe) .* ...
              (Cfe + Afe*exp(1j*phase_mod));

        % DC de ADC y ruido complejo.
        raw = raw + (0.012 - 1j*0.008);
        raw = raw + 0.003/sqrt(2) * ...
            (randn(Nfe,1) + 1j*randn(Nfe,1));

        [bb_fe, Fs_fe] = frontend_dsp(raw, cfg);

        fprintf('Fs de salida: %.3f Hz (esperado: 100 Hz)\n', Fs_fe);
        fprintf('Duracion salida: %.3f s\n', numel(bb_fe)/Fs_fe);

        % Diagnostico simple de la modulacion recuperada.
        % Se usa el centro conocido solo para esta prueba del frontend.
        % No forma parte del algoritmo real.
        centro_est = fit_circle_kasa_test(bb_fe);
        zfe = bb_fe - centro_est;
        phi_fe = unwrap(angle(zfe));
        phi_fe = detrend(phi_fe, 1);

        [fpk_fe, ~] = pico_en_banda(phi_fe, Fs_fe, 0.5, 2.0);
        fprintf('Modulacion recuperada: %.3f Hz (verdad %.3f Hz)\n\n', ...
            fpk_fe, f_mod);

        figure('Name','Test sintetico - Frontend');
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

        nexttile;
        nshow = min(numel(raw), round(0.002*cfg.fs));
        plot((0:nshow-1)/cfg.fs*1e3, real(raw(1:nshow)));
        grid on;
        xlabel('ms'); ylabel('I');
        title('RX Low-IF sintetica');

        nexttile;
        plot(real(bb_fe), imag(bb_fe), '.');
        axis equal; grid on;
        xlabel('I'); ylabel('Q');
        title('IQ tras frontend');

        nexttile;
        tbb = (0:numel(phi_fe)-1)/Fs_fe;
        plot(tbb, phi_fe);
        grid on;
        xlabel('s'); ylabel('rad');
        title('Fase recuperada');

        nexttile;
        [fax, P] = espectro_local(phi_fe, Fs_fe);
        plot(fax, 20*log10(P + 1e-12));
        xlim([0 3]); grid on;
        xline(f_mod,'--','verdad');
        xlabel('Hz'); ylabel('dB');
        title(sprintf('Pico %.3f Hz', fpk_fe));
        drawnow;

        clear raw tfe phase_mod;
    end
end

%% =========================================================================
% B) PRUEBAS DEL PROCESADOR DE SIGNOS VITALES
% =========================================================================
if exist('procesar_signos_vitales_fixed', 'file') ~= 2
    error(['No se encontro procesar_signos_vitales_fixed.m. ' ...
        'El nombre del archivo debe coincidir exactamente con la funcion.']);
end

fprintf('--- B) PRUEBAS DE RPM / BPM ---\n');

% Parametros comunes del modelo.
base = struct();
base.Fs = Fs_bb;
base.duracion_s = duracion_s;
base.fc = fc;
base.resp_mm = 4.0;              % amplitud pico de desplazamiento
base.hr_mm = 0.30;               % amplitud pico de desplazamiento
base.echo_amp = 0.030;           % radio aproximado del arco IQ
base.center_ratio = 12;          % |C| / A
base.noise_sigma = 0.0015;
base.lo_offset_hz = 0.018;
base.lo_drift_rad = 0.20;        % deriva lenta pico
base.motion_rad = 0;
base.motion_time_s = 15;
base.motion_width_s = 0.30;

casos = repmat(base, 6, 1);

casos(1).nombre = 'Nominal';
casos(1).rpm = 15;
casos(1).bpm = 72;

casos(2).nombre = 'Respiracion lenta';
casos(2).rpm = 6;
casos(2).bpm = 72;

% Caso deliberadamente dificil:
% 24 rpm = 0.4 Hz; 3er armonico = 1.2 Hz = 72 bpm.
casos(3).nombre = 'HR coincide con 3xRR';
casos(3).rpm = 24;
casos(3).bpm = 72;

casos(4).nombre = 'Bajo SNR / clutter alto';
casos(4).rpm = 15;
casos(4).bpm = 90;
casos(4).center_ratio = 45;
casos(4).noise_sigma = 0.006;
casos(4).hr_mm = 0.20;

casos(5).nombre = 'Movimiento puntual';
casos(5).rpm = 15;
casos(5).bpm = 72;
casos(5).motion_rad = 2.5;

% Control negativo: no hay movimiento respiratorio ni cardiaco.
casos(6).nombre = 'Sin signos vitales';
casos(6).rpm = NaN;
casos(6).bpm = NaN;
casos(6).resp_mm = 0;
casos(6).hr_mm = 0;
casos(6).noise_sigma = 0.002;

ncasos = numel(casos);

nombre = strings(ncasos,1);
rpm_true = nan(ncasos,1);
rpm_est = nan(ncasos,1);
err_rpm = nan(ncasos,1);
bpm_true = nan(ncasos,1);
bpm_est = nan(ncasos,1);
err_bpm = nan(ncasos,1);
rr_fallback = false(ncasos,1);
rr_low_conf = false(ncasos,1);
rr_candidate = nan(ncasos,1);
hr_fallback = false(ncasos,1);
hr_low_conf = false(ncasos,1);
hr_candidate = nan(ncasos,1);
center_ratio_out = nan(ncasos,1);
arc_span = nan(ncasos,1);
hr_agreement = nan(ncasos,1);

datos = cell(ncasos,1);
infos = cell(ncasos,1);

for k = 1:ncasos
    p = casos(k);

    fprintf('\n============================================================\n');
    fprintf('CASO %d/%d: %s\n', k, ncasos, p.nombre);
    if isfinite(p.rpm)
        fprintf('Verdad: RR %.1f rpm | HR %.1f bpm\n', p.rpm, p.bpm);
    else
        fprintf('Verdad: SIN respiracion / SIN latido\n');
    end
    fprintf('============================================================\n');

    [iq, verdad] = generar_iq_sintetico(p);

    [rr_hat, hr_hat, info] = procesar_signos_vitales_fixed( ...
        iq, Fs_bb, dHP, dLP_rr);

    nombre(k) = string(p.nombre);
    rpm_true(k) = p.rpm;
    bpm_true(k) = p.bpm;
    rpm_est(k) = rr_hat;
    bpm_est(k) = hr_hat;

    if isfinite(p.rpm) && isfinite(rr_hat)
        err_rpm(k) = rr_hat - p.rpm;
    end
    if isfinite(p.bpm) && isfinite(hr_hat)
        err_bpm(k) = hr_hat - p.bpm;
    end

    if isfield(info,'rr_fallback'), rr_fallback(k) = logical(info.rr_fallback); end
    if isfield(info,'rr_low_confidence'), rr_low_conf(k) = logical(info.rr_low_confidence); end
    if isfield(info,'rr_candidate_rpm'), rr_candidate(k) = info.rr_candidate_rpm; end
    if isfield(info,'hr_fallback'), hr_fallback(k) = logical(info.hr_fallback); end
    if isfield(info,'hr_low_confidence'), hr_low_conf(k) = logical(info.hr_low_confidence); end
    if isfield(info,'hr_candidate_bpm'), hr_candidate(k) = info.hr_candidate_bpm; end
    if isfield(info,'center_offset_ratio'), center_ratio_out(k) = info.center_offset_ratio; end
    if isfield(info,'arc_span_rad'), arc_span(k) = info.arc_span_rad; end
    if isfield(info,'hr_agreement'), hr_agreement(k) = info.hr_agreement; end

    datos{k} = struct('iq',iq,'verdad',verdad);
    infos{k} = info;

    fprintf('\nRESULTADO DEL CASO:\n');
    if isfinite(p.rpm)
        fprintf('  RR: verdad %.1f | estimado %.2f | error %+.2f rpm\n', ...
            p.rpm, rr_hat, err_rpm(k));
        fprintf('  HR: verdad %.1f | estimado %.2f | error %+.2f bpm\n', ...
            p.bpm, hr_hat, err_bpm(k));
    else
        fprintf('  RR estimado: %.2f | candidato %.2f | lowconf=%d | fallback=%d\n', ...
            rr_hat, rr_candidate(k), rr_low_conf(k), rr_fallback(k));
        fprintf('  HR estimado: %.2f | candidato %.2f | lowconf=%d | fallback=%d\n', ...
            hr_hat, hr_candidate(k), hr_low_conf(k), hr_fallback(k));
    end
end

pass_rr = false(ncasos,1);
pass_hr = false(ncasos,1);
for k = 1:ncasos
    if isfinite(rpm_true(k))
        pass_rr(k) = isfinite(rpm_est(k)) && abs(rpm_est(k)-rpm_true(k)) <= 1.0;
    else
        pass_rr(k) = ~isfinite(rpm_est(k));
    end

    if isfinite(bpm_true(k))
        % En el caso de movimiento se acepta rechazo por baja confianza.
        if k == 5
            pass_hr(k) = (~isfinite(bpm_est(k)) && hr_low_conf(k)) || ...
                         (isfinite(bpm_est(k)) && abs(bpm_est(k)-bpm_true(k)) <= 5.0);
        else
            pass_hr(k) = isfinite(bpm_est(k)) && abs(bpm_est(k)-bpm_true(k)) <= 5.0;
        end
    else
        pass_hr(k) = ~isfinite(bpm_est(k));
    end
end

T = table(nombre, rpm_true, rpm_est, err_rpm, rr_candidate, ...
    bpm_true, bpm_est, err_bpm, hr_candidate, ...
    rr_fallback, rr_low_conf, hr_fallback, hr_low_conf, ...
    center_ratio_out, arc_span, hr_agreement, pass_rr, pass_hr);

fprintf('\n\n================ RESUMEN =================\n');
disp(T);
fprintf('Casos RR aprobados: %d/%d\n', nnz(pass_rr), ncasos);
fprintf('Casos HR aprobados: %d/%d\n', nnz(pass_hr), ncasos);
if all(pass_rr) && all(pass_hr)
    fprintf('RESULTADO GLOBAL: PASS\n');
else
    fprintf('RESULTADO GLOBAL: REVISAR CASOS FALLIDOS\n');
end

%% -------------------------------------------------------------------------
% GRAFICA DE UN CASO
% -------------------------------------------------------------------------
caso_grafica = 1;  % cambiar 1..6

iq = datos{caso_grafica}.iq;
verdad = datos{caso_grafica}.verdad;
info = infos{caso_grafica};

figure('Name', sprintf('Caso sintetico: %s', casos(caso_grafica).nombre), ...
       'Units','normalized','Position',[0.07 0.08 0.86 0.80]);

tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile;
plot(real(iq), imag(iq), '.','MarkerSize',3);
axis equal; grid on;
xlabel('I'); ylabel('Q');
title('IQ sintetico');

nexttile;
plot(verdad.t, verdad.x_resp_mm, 'LineWidth',1);
hold on;
plot(verdad.t, verdad.x_hr_mm, 'LineWidth',1);
grid on;
xlabel('s'); ylabel('mm');
legend('respiracion','corazon','Location','best');
title('Desplazamientos verdaderos');

nexttile;
plot(verdad.t, verdad.phi_total, 'LineWidth',1);
grid on;
xlabel('s'); ylabel('rad');
title('Fase fisiologica verdadera');

nexttile;
if isfield(info,'x_rr') && ~isempty(info.x_rr)
    tx = (0:numel(info.x_rr)-1)/Fs_bb;
    plot(tx, info.x_rr);
    grid on;
    xlabel('s'); ylabel('rad');
    title(sprintf('RR procesada: %.2f rpm', rpm_est(caso_grafica)));
else
    axis off;
    text(0.5,0.5,'x_rr no disponible','HorizontalAlignment','center');
end

nexttile;
if isfield(info,'f_rr_axis') && ~isempty(info.f_rr_axis)
    plot(info.f_rr_axis, 20*log10(info.P_rr + 1e-12));
    xlim([0 1]); grid on;
    if isfinite(casos(caso_grafica).rpm)
        xline(casos(caso_grafica).rpm/60,'--','verdad');
    end
    xlabel('Hz'); ylabel('dB');
    title('Espectro RR');
else
    axis off;
end

nexttile;
if isfield(info,'f_hr_axis') && ~isempty(info.f_hr_axis)
    plot(info.f_hr_axis, 20*log10(info.P_hr + 1e-12));
    xlim([0.5 3]); grid on;
    if isfinite(casos(caso_grafica).bpm)
        xline(casos(caso_grafica).bpm/60,'--','verdad');
    end
    xlabel('Hz'); ylabel('dB');
    title('Espectro HR');
else
    axis off;
end

sgtitle(tl, sprintf('%s | RR %.2f rpm | HR %.2f bpm', ...
    casos(caso_grafica).nombre, rpm_est(caso_grafica), bpm_est(caso_grafica)), ...
    'FontWeight','bold');

%% Guardar datos para poder repetir procesamiento sin regenerar ruido
save('banco_sintetico_radar.mat', 'casos', 'datos', 'T', 'Fs_bb', 'fc');
fprintf('\nSe guardo banco_sintetico_radar.mat\n');

%% =========================================================================
% FUNCIONES LOCALES
% =========================================================================
function [iq, verdad] = generar_iq_sintetico(p)

Fs = p.Fs;
N = round(p.duracion_s * Fs);
t = (0:N-1).' / Fs;

c = 299792458;
lambda = c / p.fc;

if isfinite(p.rpm)
    f_rr = p.rpm / 60;
else
    f_rr = 0;
end
if isfinite(p.bpm)
    f_hr = p.bpm / 60;
else
    f_hr = 0;
end

% Respiracion: fundamental + pequena segunda componente para evitar una
% senal demasiado ideal.
if f_rr > 0 && p.resp_mm > 0
    x_resp = p.resp_mm*1e-3 * ...
        (sin(2*pi*f_rr*t) + 0.10*sin(4*pi*f_rr*t + 0.4));
else
    x_resp = zeros(N,1);
end

% Latido: desplazamiento mucho menor y ligeramente no sinusoidal.
if f_hr > 0 && p.hr_mm > 0
    x_hr = p.hr_mm*1e-3 * ...
        (sin(2*pi*f_hr*t + 0.7) + ...
         0.22*sin(4*pi*f_hr*t + 1.1));
else
    x_hr = zeros(N,1);
end

% Movimiento corporal puntual. Se expresa como fase adicional porque un
% movimiento grande no necesita parecerse a una oscilacion fisiologica.
if p.motion_rad ~= 0
    phi_motion = p.motion_rad * ...
        exp(-0.5*((t - p.motion_time_s)/p.motion_width_s).^2);
else
    phi_motion = zeros(N,1);
end

phi_resp = 4*pi*x_resp/lambda;
phi_hr   = 4*pi*x_hr/lambda;
phi_phys = phi_resp + phi_hr + phi_motion;

% Centro IQ: fuga TX->RX + clutter.
A = p.echo_amp;
C = p.center_ratio * A * exp(1j*0.35);

% Offset constante + deriva lenta del LO.
phi_lo = 2*pi*p.lo_offset_hz*t;
if p.lo_drift_rad ~= 0
    phi_lo = phi_lo + p.lo_drift_rad * sin(2*pi*0.015*t + 0.2);
end

iq_clean = exp(1j*phi_lo) .* ...
    (C + A*exp(1j*(phi_phys + 0.8)));

noise = p.noise_sigma/sqrt(2) * ...
    (randn(N,1) + 1j*randn(N,1));

iq = iq_clean + noise;

verdad = struct();
verdad.t = t;
verdad.lambda = lambda;
verdad.x_resp_mm = x_resp*1e3;
verdad.x_hr_mm = x_hr*1e3;
verdad.phi_resp = phi_resp;
verdad.phi_hr = phi_hr;
verdad.phi_motion = phi_motion;
verdad.phi_total = phi_phys;
verdad.phi_lo = phi_lo;
verdad.iq_clean = iq_clean;
end


function c = fit_circle_kasa_test(z)
% Centro algebraico solo para diagnosticar la prueba del frontend.
x = real(z(:));
y = imag(z(:));
A = [x, y, ones(numel(x),1)];
b = x.^2 + y.^2;
coef = A \ b;
c = complex(coef(1)/2, coef(2)/2);
end

function [fpeak, peakval] = pico_en_banda(x, Fs, fmin, fmax)
[f, P] = espectro_local(x, Fs);
m = f >= fmin & f <= fmax;
if ~any(m)
    fpeak = NaN;
    peakval = NaN;
    return;
end
idx = find(m);
[peakval, ir] = max(P(idx));
fpeak = f(idx(ir));
end

function [f, P] = espectro_local(x, Fs)
x = x(:);
N = numel(x);
Nfft = 2^nextpow2(max(4096, 8*N));
w = hamming(N);
X = fft((x-mean(x)).*w, Nfft);
K = floor(Nfft/2)+1;
P = abs(X(1:K));
f = (0:K-1).' * Fs/Nfft;
end
