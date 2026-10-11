% TEST_SINTETICO Prueba del procesador sin bladeRF.
clear; clc; close all;
rng(1);

Fs = 100;
T = 30;
t = (0:1/Fs:T-1/Fs).';

rpm_true = 15;
bpm_true = 72;
f_rr = rpm_true/60;
f_hr = bpm_true/60;

% Movimiento equivalente expresado directamente como fase radar.
phi_resp = 1.00*sin(2*pi*f_rr*t);
phi_heart = 0.10*sin(2*pi*f_hr*t + 0.4);
phi = phi_resp + phi_heart;

% Arco IQ con clutter estatico y pequeno offset residual TX/RX.
c = 0.60 + 0.25j;
R = 0.12;
f_off = 0.020;
iq = (c + R*exp(1j*phi)) .* exp(1j*2*pi*f_off*t);
iq = iq + 0.002*(randn(size(iq)) + 1j*randn(size(iq)));

dHP = designfilt('highpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',0.04,'SampleRate',Fs);
dLP_rr = designfilt('lowpassiir','FilterOrder',4, ...
    'HalfPowerFrequency',0.8,'SampleRate',Fs);

[rpm,bpm,info] = procesar_signos_vitales_fixed(iq,Fs,dHP,dLP_rr);

fprintf('\nVerdad     : RR %.1f rpm | HR %.1f bpm\n',rpm_true,bpm_true);
fprintf('Estimacion : RR %.1f rpm | HR %.1f bpm\n',rpm,bpm);
fprintf('Confianza  : RR low=%d | HR low=%d\n', ...
    info.rr_low_confidence,info.hr_low_confidence);

pass_rr = isfinite(rpm) && abs(rpm-rpm_true) <= 1.0;
pass_hr = isfinite(bpm) && abs(bpm-bpm_true) <= 5.0;
if pass_rr && pass_hr
    fprintf('TEST: PASS\n');
else
    warning('TEST: FAIL. Revise la cadena antes de probar con hardware.');
end
