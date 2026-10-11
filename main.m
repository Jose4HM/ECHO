% RADAR CW DE SIGNOS VITALES - SCRIPT PRINCIPAL
clear; clc; close all;

% -------------------------------------------------------------------------
% 1. CAPTURA
% -------------------------------------------------------------------------
try
    [rx_continuous, cfg] = cw();
catch ME
    error('Fallo en cw.m: %s', ME.message);
end

if isempty(rx_continuous)
    error('cw.m no genero muestras.');
end
if isfield(cfg, 'overruns') && cfg.overruns > 0
    error(['La captura tuvo %d overruns. No se procesa porque faltan muestras ' ...
        'y la fase deja de ser continua.'], cfg.overruns);
end

% -------------------------------------------------------------------------
% 2. LOW-IF -> BANDA BASE
% -------------------------------------------------------------------------
fprintf('\nEjecutando front-end DSP...\n');
[iq_frontend, Fs_frontend] = frontend_dsp(rx_continuous, cfg);
if isempty(iq_frontend) || ~isfinite(Fs_frontend) || Fs_frontend <= 0
    error('frontend_dsp devolvio una salida invalida.');
end

% -------------------------------------------------------------------------
% 3. FILTROS FISIOLOGICOS
% -------------------------------------------------------------------------
% Un unico paso-alto para todo el rango respiratorio de interes.
% fc_hp = 0.04 Hz queda suficientemente por debajo del limite inferior
% del detector de RR (6 rpm = 0.10 Hz), evitando atenuar respiraciones lentas
% sin necesidad de seleccionar previamente un modo de respiracion.
fc_hp = 0.04;

dur_s = numel(iq_frontend)/Fs_frontend;
if dur_s < 25
    warning(['Captura de solo %.1f s con HP %.2f Hz. Para estimaciones mas ' ...
        'estables de RR/HR se recomiendan aproximadamente 30 s.'], dur_s, fc_hp);
end

fprintf('Filtro fisiologico: HP %.2f Hz, duracion %.1f s.\n', fc_hp, dur_s);

dHP = designfilt('highpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', fc_hp, 'SampleRate', Fs_frontend);

dLP_rr = designfilt('lowpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency', 0.8, 'SampleRate', Fs_frontend);

% -------------------------------------------------------------------------
% 4. ESTIMACION
% -------------------------------------------------------------------------
[rpm, bpm, vital_info] = procesar_signos_vitales_fixed( ...
    iq_frontend, Fs_frontend, dHP, dLP_rr);

fprintf('\n=============================\n');
fprintf('RESULTADOS\n');
if vital_info.rr_fallback || ~isfinite(rpm)
    fprintf('Respiracion: SIN DETECCION\n');
elseif vital_info.rr_low_confidence
    fprintf('Respiracion: %.1f rpm  [BAJA CONFIANZA]\n', rpm);
else
    fprintf('Respiracion: %.1f rpm\n', rpm);
end

if vital_info.hr_fallback || ~isfinite(bpm)
    fprintf('Ritmo cardiaco: SIN DETECCION\n');
elseif vital_info.hr_low_confidence
    fprintf('Ritmo cardiaco: %.1f bpm  [BAJA CONFIANZA]\n', bpm);
else
    fprintf('Ritmo cardiaco: %.1f bpm\n', bpm);
end

if isfinite(vital_info.hr_agreement)
    fprintf('Acuerdo HR entre ventanas: %.0f%%\n', 100*vital_info.hr_agreement);
end
fprintf('=============================\n');

plot_diagnostico_radar(rx_continuous, cfg, iq_frontend, Fs_frontend, ...
    rpm, bpm, vital_info);

% =========================================================================
function plot_diagnostico_radar(rx_raw, cfg, iq_bb, Fs_bb, rpm, bpm, info)
Fs_raw = cfg.fs;
rx_raw = rx_raw(:);
iq_bb = iq_bb(:);

rr_txt = 'SIN DET.';
if isfinite(rpm), rr_txt = sprintf('%.1f rpm', rpm); end
hr_txt = 'SIN DET.';
if isfinite(bpm), hr_txt = sprintf('%.1f bpm', bpm); end

figure('Name', 'Diagnostico radar signos vitales', ...
    'Units', 'normalized', 'Position', [0.04 0.06 0.92 0.86]);
sgtitle(sprintf('Radar | RR %s | HR %s | |c|/R %.0f | arco %.2f rad', ...
    rr_txt, hr_txt, info.center_offset_ratio, info.arc_span_rad), ...
    'FontWeight', 'bold', 'FontSize', 13);

% 1) Nivel RX por bloque
subplot(3,3,1);
if isfield(cfg,'block_peak') && isfield(cfg,'block_rms')
    n_blk = numel(cfg.block_peak);
    t_blk = (0:n_blk-1)*(cfg.block_len/cfg.fs);
    plot(t_blk,cfg.block_peak,'b'); hold on;
    plot(t_blk,cfg.block_rms,'c');
    yline(0.99,':r');
    idx_inc = cfg.actual_counts ~= cfg.block_len;
    if any(idx_inc)
        plot(t_blk(idx_inc),cfg.block_peak(idx_inc),'rx');
    end
    grid on; xlabel('s'); ylabel('amplitud');
    title(sprintf('RX (%d bloques incompletos)', nnz(idx_inc)));
else
    n_raw = min(numel(rx_raw), round(0.01*Fs_raw));
    t_ms = (0:n_raw-1)/Fs_raw*1e3;
    plot(t_ms,real(rx_raw(1:n_raw))); hold on;
    plot(t_ms,imag(rx_raw(1:n_raw)));
    grid on; xlabel('ms'); ylabel('amplitud'); title('RX cruda');
end

% 2) IQ sin centrar
subplot(3,3,2);
idx_c = round(linspace(1,numel(iq_bb),min(8000,numel(iq_bb))));
plot(real(iq_bb(idx_c)),imag(iq_bb(idx_c)),'.','MarkerSize',3); hold on;
plot(0,0,'r+','MarkerSize',10,'LineWidth',1.5);
axis equal; grid on; xlabel('I'); ylabel('Q');
title(sprintf('IQ banda base @ %.0f Hz', Fs_bb));

% 3) FFT cruda low-IF
subplot(3,3,3);
plot_fft_complex(rx_raw,Fs_raw,262144,[-150 150],1e3);
xline(cfg.tx_tone_hz/1e3,'--r',sprintf('%.0f kHz',cfg.tx_tone_hz/1e3));
grid on; xlabel('kHz'); ylabel('dB'); title('FFT RX cruda');

% 4) Fase
subplot(3,3,4);
t = (0:numel(info.phi)-1).'/info.Fs;
plot(t,info.phi,'LineWidth',0.8);
grid on; xlabel('s'); ylabel('rad'); title('Fase desenvuelta');

% 5) Respiracion
subplot(3,3,5);
t = (0:numel(info.x_rr)-1).'/info.Fs;
plot(t,info.x_rr,'LineWidth',1.0);
grid on; xlabel('s'); ylabel('rad'); title(['Respiracion: ' rr_txt]);

% 6) Cardiaca
subplot(3,3,6);
t = (0:numel(info.x_hr)-1).'/info.Fs;
plot(t,info.x_hr,'LineWidth',0.8);
grid on; xlabel('s'); ylabel('rad'); title(['Cardiaca: ' hr_txt]);

% 7) Espectro RR
subplot(3,3,7);
plot(info.f_rr_axis,20*log10(info.P_rr+1e-12),'LineWidth',1.1); hold on;
if isfinite(info.f_rr_est)
    xline(info.f_rr_est,'--r',rr_txt);
end
grid on; xlim([0 1.2]); xlabel('Hz'); ylabel('dB'); title('Espectro respiratorio');

% 8) Espectro HR
subplot(3,3,8);
plot(info.f_hr_axis,20*log10(info.P_hr+1e-12),'LineWidth',1.1); hold on;
if isfinite(info.f_hr_est)
    xline(info.f_hr_est,'--r',hr_txt);
end
if isfinite(info.f_rr_est)
    for h = 2:6
        fh = h*info.f_rr_est;
        if fh >= 0.4 && fh <= 3.6
            xline(fh,':',sprintf('%dxRR',h),'Color',[0.6 0.6 0.6]);
        end
    end
end
grid on; xlim([0.4 3.6]); xlabel('Hz'); ylabel('dB'); title('Espectro cardiaco');

% 9) HR por ventanas
subplot(3,3,9);
if ~isempty(info.hr_window_bpm)
    stem(1:numel(info.hr_window_bpm),info.hr_window_bpm,'filled'); hold on;
    if isfinite(bpm), yline(bpm,'--r'); end
    grid on; xlabel('ventana'); ylabel('bpm');
    if isfinite(info.hr_agreement)
        title(sprintf('HR por ventana (acuerdo %.0f%%)',100*info.hr_agreement));
    else
        title('HR por ventana');
    end
else
    text(0.5,0.5,'sin ventanas','HorizontalAlignment','center'); axis off;
end
end

function plot_fft_complex(x, Fs, nfft_max, xlim_values, unit_scale)
x = x(:);
N = min(numel(x),nfft_max);
Nfft = 2^nextpow2(N);
w = hamming(N);
X = fftshift(fft(x(1:N).*w,Nfft));
f = (-Nfft/2:Nfft/2-1)*Fs/Nfft;
plot(f/unit_scale,20*log10(abs(X)+1e-12),'LineWidth',1.0);
xlim(xlim_values);
end
