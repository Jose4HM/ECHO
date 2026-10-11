function [rpm_est, bpm_est, info] = procesar_signos_vitales_fixed(buffer_IQ, Fs, dHP, dLP_rr, rango_s)
%PROCESAR_SIGNOS_VITALES_FIXED  Respiracion y ritmo cardiaco desde IQ de un
%   radar CW Doppler (bladeRF 2.0), en banda base a ~100 Hz.
%
% Entradas:
%   buffer_IQ – Vector complejo (Nx1) a Fs, CON su offset de banda base
%               intacto (fuga TX->RX + clutter). NO se lo quites antes.
%   Fs        – Frecuencia de muestreo del buffer (Hz). Se espera ~100 Hz.
%   dHP       – Paso-alto IIR ya disenado a Fs (quita deriva lenta)
%   dLP_rr    – Paso-bajo IIR ya disenado a Fs (banda de respiracion)
%
% Salidas:
%   rpm_est – Tasa respiratoria estimada (rpm)
%   bpm_est – Frecuencia cardiaca estimada (bpm)
%   info    – Struct de diagnostico. Campos clave:
%       .center_offset_ratio  |c|/R : cuanto pesa la fuga+clutter estatico
%                             frente a la senal util. <20 excelente,
%                             20-100 aceptable, >100 hay que mejorar el
%                             aislamiento o la geometria.
%       .arc_span_rad         Apertura angular del arco IQ (rad). La
%                             respiracion deberia dar 0.5-3 rad.
%       .hr_window_bpm        HR estimado en cada ventana deslizante
%       .hr_agreement         Fraccion de ventanas que coinciden +-5 bpm
%       .rr_fallback / .hr_fallback / .hr_low_confidence


% ---------------- Parametros ---------------------------------------------
pad_mult      = 8;     % Zero-padding para interpolacion espectral fina
% Segundos a descartar en cada borde tras filtrar. Se adapta a la duracion:
% en un tramo corto (apnea) recortar 5 s por lado se lleva la mitad de los
% datos. Nunca menos de 1.5 s, que es el transitorio real de filtfilt.
trim_s        = 5;
circ_win_s    = 15;    % Ventana para el ajuste de circulo
sic_win_s     = 10;    % Ventana para la cancelacion de respiracion
sic_nharm     = 6;     % Armonicos de respiracion a cancelar
hr_win_s      = 20;    % Ventana para el chequeo de consistencia de HR
welch_win_s   = 15;    % Ventana del espectro promediado para HR

% ---------------- Adaptacion a la duracion -------------------------------
% Todos los valores de arriba estaban dimensionados para capturas de 60 s.
% La referencia (VitRad, HardwareX 2022) graba 8 s, y una captura corta es
% MEJOR aqui: la deriva entre los dos sintetizadores del AD9361 se acumula
% con el tiempo (medimos ~0.06 vueltas de fase por segundo), asi que 15 s
% dejan ~0.9 vueltas frente a las 3.5 de una captura de 60 s. Tambien baja
% la probabilidad de que el sujeto se mueva.
%
% Con menos datos hay que reescalar las ventanas o los estimadores se
% quedan con una sola y dejan de promediar.
dur_cap_s = size(buffer_IQ(:),1) / Fs;

if dur_cap_s < 45
    trim_s      = max(1.2, 0.10 * dur_cap_s);
    circ_win_s  = max(5.0, 0.45 * dur_cap_s);
    % La ventana de Welch NO debe escalar con la captura: la fija la
    % variabilidad del ritmo cardiaco, que mantiene la frecuencia estable
    % unos 10-15 s. Mas larga smearea el pico; mas corta pierde
    % resolucion. Se acota a esa franja para que una captura de 30 s
    % rinda 5-6 ventanas promediadas en vez de 2.
    welch_win_s = min(14.0, max(7.0, 0.40 * dur_cap_s));
    hr_win_s    = max(5.0, 0.40 * dur_cap_s);
    sic_win_s   = max(5.0, 0.35 * dur_cap_s);
    fprintf(['Captura de %.0f s: ventanas reescaladas ' ...
        '(circulo %.1f s, Welch %.1f s, consistencia %.1f s).\n'], ...
        dur_cap_s, circ_win_s, welch_win_s, hr_win_s);
end

% El cancelador de respiracion necesita varias ventanas para adaptarse. Con
% pocos datos ajusta 6 armonicos sobre una sola ventana, que es justo la
% situacion en la que se come el latido. La cadena de la referencia (su
% figura 8) NO lleva cancelador: arctan, paso-banda, ventana, zero-padding
% y FFT. Por debajo de 20 s se usa esa cadena simple.
% ####################  INTERRUPTOR DEL SIC  ####################
% Cancelacion armonica de la respiracion. NO es lo mismo que 'protocolo'
% de main.m: aquel fija el corte del paso-alto, este decide si se resta el
% modelo armonico de la respiracion antes de buscar el latido.
%
%   'auto'  se activa solo si la captura llega a 20 s
%   'si'    forzado ON
%   'no'    forzado OFF  <- la cadena de VitRad (su figura 8) NO lleva SIC
%
% Protocolo recomendado: prueba primero en 'no'. Si el latido sale, tu
% cadena coincide con la de la referencia y es mas facil de defender. Solo
% si no sale, repite el MISMO .mat en 'si' y compara: eso te da una figura
% de "con y sin cancelacion armonica" sobre datos identicos.
%
% OJO: el SIC ajusta 6 armonicos. El refinamiento de frecuencia se bajo de
% +-12% a +-3% precisamente porque con el margen ancho un armonico podia
% desplazarse hasta el latido y restarlo como si fuera respiracion. Aun
% asi, comprueba siempre que el HR resultante no caiga sobre un multiplo
% de RR.
sic_modo = 'no';        % <-- 'auto' | 'si' | 'no'

switch sic_modo
    case 'auto', usar_cancelador = dur_cap_s >= 20;
    case 'si',   usar_cancelador = true;
    case 'no',   usar_cancelador = false;
    otherwise,   error('sic_modo debe ser ''auto'', ''si'' o ''no''.');
end
if usar_cancelador
    fprintf('SIC: modo ''%s'' -> cancelador ACTIVADO.\n', sic_modo);
else
    fprintf('SIC: modo ''%s'' -> cancelador desactivado.\n', sic_modo);
end
% ##############################################################

info = struct('rr_fallback', false, 'hr_fallback', false, ...
    'hr_low_confidence', false, 'hr_prominence', 0, 'hr_score_ratio', 0, ...
    'center_offset_ratio', NaN, 'arc_span_rad', NaN, 'circle_radius', NaN, ...
    'hr_window_bpm', [], 'hr_agreement', NaN, 'Fs', NaN, ...
    'phi', [], 'x', [], 'x_rr', [], 'x_hr', [], ...
    'f_rr_axis', [], 'P_rr', [], 'f_hr_axis', [], 'P_hr', [], ...
    'f_rr_est', NaN, 'f_hr_est', NaN);

if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('Fs debe ser un escalar positivo en Hz.');
end
info.Fs = Fs;

if Fs > 500
    warning(['Fs = %.1f Hz es innecesariamente alta para senales de <4 Hz. ' ...
        'Los filtros fisiologicos quedan mal condicionados. Usa el ' ...
        'frontend v4, que entrega 100 Hz.'], Fs);
end

buffer_IQ = buffer_IQ(:);

check_filter_fs(dHP, 'dHP', Fs);
check_filter_fs(dLP_rr, 'dLP_rr', Fs);

% ================================================================
% 1. ACONDICIONAMIENTO IQ
% ================================================================
if any(~isfinite(buffer_IQ))
    n_bad = sum(~isfinite(buffer_IQ));
    warning('buffer_IQ contiene %d muestras no finitas (NaN/Inf); se interpolaran.', n_bad);
    idx_good = isfinite(buffer_IQ);
    if sum(idx_good) < 2
        error('buffer_IQ no tiene suficientes muestras finitas para procesar.');
    end
    idx_all = (1:numel(buffer_IQ)).';
    buffer_IQ = interp1(idx_all(idx_good), buffer_IQ(idx_good), idx_all, 'linear', 'extrap');
end

N = numel(buffer_IQ);
if N < 20 * Fs
    warning('buffer_IQ tiene %.2f s de datos; BPM/RPM seran poco estables. Usa >=40 s.', N / Fs);
end

% --- ACONDICIONAMIENTO I/Q -----------------------------------------------
% Gram-Schmidt DESACTIVADO.
S_iq = buffer_IQ;
% ================================================================
% 1b. PORTADORA RESIDUAL (offset entre el LO de TX y el de RX)
% ================================================================
% En el AD9361 los sintetizadores de TX y RX son PLL fraccionales
% independientes. Aunque compartan referencia, el redondeo fraccional deja
% un offset de fracciones de Hz, que en 50 s se acumula en varias vueltas
% completas de fase.
%
% Esto TIENE que quitarse ANTES del ajuste de circulo. Si no, la
% trayectoria IQ es un circulo grande centrado en el ORIGEN (radio =
% |clutter|) con la senal del pecho como un epiciclo minusculo encima.
% Kasa ajusta ese circulo, devuelve centro ~= 0, y angle() acaba midiendo
% la rotacion del LO con el pecho atenuado por A_pecho/A_clutter. El
% detrend polinomico posterior quita la rampa pero NO recupera esa
% atenuacion: el dano ya esta hecho.
% Con capturas cortas se omite el seguimiento de deriva (etapa 3): un
% filtro de 0.06 Hz sobre 15 s queda dominado por sus propios bordes, y con
% tan poca deriva acumulada el detrend polinomico posterior ya la absorbe.
if dur_cap_s >= 30
    [S_iq, f_off, vueltas, deriva_pp] = quitar_offset_portadora(S_iq, Fs);
else
    [S_iq, f_off, vueltas, deriva_pp] = quitar_offset_portadora(S_iq, Fs, 0);
end
info.carrier_offset_hz = f_off;

fprintf('\n--- PORTADORA RESIDUAL ---\n');
fprintf('Offset LO TX-RX estimado    : %+.4f Hz\n', f_off);
fprintf('Vueltas de fase acumuladas  : %.2f en %.0f s\n', ...
    vueltas, numel(S_iq)/Fs);
fprintf('Deriva lenta retirada       : %.2f rad pico a pico\n', deriva_pp);
info.lo_drift_pp_rad = deriva_pp;
if deriva_pp > 25
    warning(['Deriva de LO de %.1f rad. Es mucho incluso para dos PLL ' ...
        'independientes: revisa que el bladeRF no este cambiando de ' ...
        'temperatura durante la captura (dejalo 5 min encendido antes).'], ...
        deriva_pp);
end
fprintf('--------------------\n');

% ================================================================
% 2. CENTRADO DEL CIRCULO Y EXTRACCION DE FASE
% ================================================================
[c_t, c0, R0, n_fits] = estimar_centro_circulo(S_iq, Fs, circ_win_s);

info.circle_radius = R0;
info.center_offset_ratio = abs(c0) / (R0 + eps);

Z = S_iq - c_t;             % <-- centrado CORRECTO (centro del circulo)

phi_raw = unwrap(angle(Z));

info.arc_span_rad = max(phi_raw) - min(phi_raw);

fprintf('\n--- GEOMETRIA IQ ---\n');
fprintf('Radio del arco R            : %.5f\n', R0);
fprintf('|centro| / R                : %.1f\n', info.center_offset_ratio);
fprintf('Apertura angular del arco   : %.2f rad (%.0f grados)\n', ...
    info.arc_span_rad, rad2deg(info.arc_span_rad));
fprintf('Ventanas usadas en el ajuste: %d\n', n_fits);
if info.center_offset_ratio > 100
    warning(['|c|/R = %.0f : la fuga TX->RX + clutter estatico domina muy por ' ...
        'encima de la senal util. Separa mas las antenas, revisa la ' ...
        'polarizacion (TX y RX del MISMO sentido circular se rechazan el ' ...
        'eco especular) o acerca al sujeto.'], info.center_offset_ratio);
end
if info.arc_span_rad < 0.25
    warning(['Apertura del arco = %.2f rad. Muy poca modulacion: el sujeto ' ...
        'apenas mueve el pecho en la linea de vision, o las antenas no ' ...
        'apuntan al esternon.'], info.arc_span_rad);
end
if info.arc_span_rad > 6.0
    warning(['Apertura del arco = %.2f rad (> 2*pi). Probablemente hay ' ...
        'movimiento corporal grande, no solo respiracion. Repite la ' ...
        'captura con el sujeto quieto y la espalda apoyada.'], info.arc_span_rad);
end
fprintf('--------------------\n');

% Guardia de radio pequeno: con el centrado correcto casi nunca deberia
% activarse. Si se activa, es senal de que el ajuste de circulo fallo.
n = (0:numel(phi_raw)-1).';
radio = abs(Z);
% OJO: antes esto usaba prctile(radio,5) como umbral, lo que por definicion
% descarta el 5% de las muestras SIEMPRE, hubiera problema o no. El aviso
% saltaba en toda captura sana. Ahora el umbral es absoluto respecto a la
% mediana del radio.
radio_min = 0.20 * median(radio);
idx_ok = isfinite(phi_raw) & radio > radio_min;
frac_malo = 1 - nnz(idx_ok)/numel(idx_ok);
if frac_malo > 0.02
    warning(['%.1f%% de las muestras pasan muy cerca del origen tras centrar. ' ...
        'Eso sugiere que el ajuste de circulo no convergio bien.'], 100*frac_malo);
end
if nnz(idx_ok) >= max(16, round(0.5 * numel(phi_raw)))
    phi_raw = interp1(n(idx_ok), phi_raw(idx_ok), n, 'linear', 'extrap');
end

% ================================================================
% 3. DERIVA DE LO Y FILTRADO
% ================================================================
% Detrend polinomico de orden 3: la deriva termica entre los sintetizadores
% de TX y RX durante ~60 s no es una rampa perfecta. Orden 3 la sigue sin
% morder la banda respiratoria (0.1-0.75 Hz).
p = polyfit(n, phi_raw, 3);
phi_det = phi_raw - polyval(p, n);

x = filtfilt(dHP, phi_det);

% Recorte de bordes: elimina el transitorio de filtfilt y los efectos del
% ajuste polinomico cerca de los extremos.
dur_s  = numel(phi_det)/Fs;
trim_s = max(1.5, min(trim_s, 0.12 * dur_s));
n_trim = round(trim_s * Fs);
if numel(x) > 4*n_trim + round(10*Fs)
    x = x(n_trim+1 : end-n_trim);
    phi_raw_t = phi_raw(n_trim+1 : end-n_trim);
else
    warning('Captura demasiado corta para recortar %.1f s en cada borde; no se recorta.', trim_s);
    phi_raw_t = phi_raw;
end
x = x(:);
t0_util = trim_s;      % instante (s, respecto a la captura) de x(1)

% La banda respiratoria se filtra sobre la senal COMPLETA y luego se
% recorta, no al reves: asi el tramo de apnea se compara contra la
% respiracion real del resto de la captura, y no hay transitorios de
% filtro dentro del tramo.
x_rr_full = filtfilt(dLP_rr, x);
rr_pp_full = prctile(x_rr_full, 97.5) - prctile(x_rr_full, 2.5);

% ----------------------------------------------------------------------
% RECORTE DE ANALISIS (apnea, o saltarse un movimiento puntual)
% ----------------------------------------------------------------------
% Se aplica AQUI, no al principio. El ajuste de circulo de Kasa necesita un
% arco amplio para localizar el centro: la respiracion barre 1-3 rad y lo
% define bien, pero el latido solo, en apnea, son 0.05-0.13 rad -- un arco
% casi recto, con el que el ajuste esta mal condicionado y el centro sale
% donde sea. Lo mismo vale para el seguimiento de deriva del LO, que con
% 17-20 s queda dominado por los bordes del filtro.
%
% Por eso la geometria y la deriva se estiman con la captura ENTERA (donde
% si hay respiracion) y solo la busqueda espectral se restringe al tramo.
if nargin >= 5 && ~isempty(rango_s)
    assert(numel(rango_s)==2 && rango_s(2) > rango_s(1), ...
        'rango_s debe ser [inicio fin] en segundos de la captura original.');
    i0 = round((rango_s(1) - t0_util) * Fs) + 1;
    i1 = round((rango_s(2) - t0_util) * Fs);
    i0 = max(1, i0);
    i1 = min(numel(x), i1);
    if i1 - i0 + 1 < round(8*Fs)
        warning(['El tramo pedido deja solo %.1f s tras el recorte de bordes ' ...
            '(la senal util empieza en t = %.1f s). No se restringe.'], ...
            max(0,(i1-i0+1))/Fs, t0_util);
    else
        x         = x(i0:i1);
        x_rr_full = x_rr_full(i0:i1);
        phi_raw_t = phi_raw_t(i0:i1);
        fprintf(['Analisis restringido a %.1f-%.1f s (geometria y deriva ' ...
            'estimadas con la captura completa).\n'], rango_s(1), rango_s(2));
    end
end

info.phi = phi_raw_t;
info.x = x;

fprintf('Senal de fase util: %.1f s (%d muestras a %.1f Hz)\n', numel(x)/Fs, numel(x), Fs);

% ================================================================
% 4. ESTIMACION DE RESPIRACION (RR)
% ================================================================
x_rr = x_rr_full;
info.x_rr = x_rr;

[f_rr_ax, P_rr, bin_hz] = espectro(x_rr, Fs, pad_mult);
info.f_rr_axis = f_rr_ax;
info.P_rr = P_rr;

rr_min_hz = 6 / 60;
rr_max_hz = 45 / 60;
mask_rr = (f_rr_ax >= rr_min_hz) & (f_rr_ax <= rr_max_hz);

f_rr_est = NaN;
if any(mask_rr)
    f_rr_est = select_rr_peak_robust(f_rr_ax, P_rr, mask_rr, bin_hz);
end

if isnan(f_rr_est) || f_rr_est < rr_min_hz || f_rr_est > rr_max_hz
    info.rr_fallback = true;
    f_rr_est = 0.25;   % 15 rpm generico
end
rpm_est = double(60 * f_rr_est(1));
info.f_rr_est = f_rr_est;

% ================================================================
% 5. CANCELACION DE RESPIRACION (adaptativa por ventanas)
% ================================================================
% DETECCION DE APNEA. El cancelador ajusta 6 armonicos con refinamiento de
% frecuencia +-12%: es un modelo MUY flexible. Si de verdad hay
% respiracion, cancela respiracion. Si NO la hay, ese mismo modelo se
% ajusta a lo unico periodico que queda -- el latido -- y lo BORRA.
%
% El criterio es la fraccion de energia que vive en la banda respiratoria
% frente al total. Respirando normal supera el 80%; en apnea cae muy por
% debajo.
% CRITERIO: amplitud de la excursion respiratoria, no fraccion de energia.
% Respirar mueve el pecho 4-12 mm = 1-3 rad a lambda 50 mm. La deriva
% residual en apnea es de decimas de radian. Se mide de percentil 2.5 a
% 97.5 para no depender de un pico aislado.
%
% (El criterio anterior, "fraccion de energia por debajo de 0.8 Hz", no
% servia: ahi vive TODO lo lento -- deriva de LO, balanceo postural, linea
% base -- asi que daba ~97% incluso sin respirar.)
rr_pp = prctile(x_rr, 97.5) - prctile(x_rr, 2.5);
info.rr_pp_rad = rr_pp;

apnea = (rr_pp < 0.40) || (rr_pp < 0.30 * rr_pp_full);
if ~usar_cancelador
    fprintf(['Captura corta (%.0f s): se usa la cadena simple de la ' ...
        'referencia, sin cancelador de respiracion.\n'], dur_cap_s);
end
info.apnea_detectada = apnea;

fprintf('Excursion respiratoria: %.2f rad en el tramo', rr_pp);
if abs(rr_pp_full - rr_pp) > 1e-9
    fprintf(' (%.2f rad en la captura completa)', rr_pp_full);
end
fprintf('\n');
if apnea || ~usar_cancelador
    if apnea
        fprintf(['-> APNEA detectada. Se OMITE la cancelacion de respiracion:\n' ...
                 '   con el sujeto sin respirar, el modelo armonico se ajustaria\n' ...
                 '   al latido y lo eliminaria.\n']);
        info.rr_fallback = true;   % el RR de este tramo no significa nada
    end
    x_rrfree = x;
else
    x_rrfree = cancelar_respiracion(x, Fs, f_rr_est, sic_nharm, sic_win_s);
end

% ================================================================
% 6. ESTIMACION DE RITMO CARDIACO (HR)
% ================================================================
margin  = max(0.04, 0.08 * f_rr_est);
if apnea
    fmin_hr = 45/60;      % sin respiracion, no hay que esquivarla
else
    fmin_hr = max(45/60, f_rr_est + margin);
end
% BANDA DE BUSQUEDA CARDIACA. Antes 45-200 bpm, que es todo el rango
% clinico posible (bradicardia severa a taquicardia de esfuerzo). Para un
% adulto EN REPOSO Y SENTADO, que es el protocolo de medida, el rango
% fisiologico es 60-100 bpm. Restringirlo es informacion previa legitima y
% es lo que hace la literatura del area.
%
% Motivo concreto: con la banda ancha, el penalizador de armonicos descarta
% la meseta 0.5-1.5 Hz (donde caen los armonicos 2o a 6o de la respiracion,
% y tambien el latido) y elige un maximo local en 2.05 Hz que esta LIBRE de
% armonicos pero tambien de senal. Medido: 122.7 bpm con verdad de 72.
hr_min_bpm = 55;
hr_max_bpm = 160;
fmin_hr = max(fmin_hr, hr_min_bpm/60);
fmax_hr = hr_max_bpm/60;
fprintf('Banda cardiaca de busqueda: %.0f-%.0f bpm (adulto en reposo).\n', ...
    60*fmin_hr, 60*fmax_hr);
if fmin_hr >= fmax_hr - 0.1
    fmin_hr = hr_min_bpm/60;
end

dHR = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', fmin_hr, ...
    'HalfPowerFrequency2', fmax_hr, ...
    'SampleRate', Fs);
x_hr = filtfilt(dHR, x_rrfree);
info.x_hr = x_hr;

% Espectro PROMEDIADO (Welch) en vez de una unica FFT larga.
%
% El corazon no es un tono estable: la variabilidad normal mueve el ritmo
% +-5 bpm a lo largo del minuto. Una FFT de 50 s tiene resolucion de
% 0.02 Hz = 1.2 bpm, asi que la energia cardiaca se reparte entre ~8 bins y
% el pico se hunde ~9 dB. Los picos de ruido son puntuales y NO se
% reparten, con lo que la FFT larga favorece al ruido frente a la senal.
%
% Con ventanas de 15 s la resolucion (4 bpm) coincide con la anchura real
% del pico cardiaco, y promediar ~10 ventanas solapadas reduce la varianza
% del ruido sin tocar la senal.
[f_hr_ax, P_hr] = espectro_welch(x_hr, Fs, welch_win_s, 0.25, pad_mult);

% ================= BLANQUEO ESPECTRAL =================
% El espectro cardiaco tiene un hombro grande en 0.8-1.4 Hz (armonicos de
% respiracion) y cae ~40 dB hacia 3.5 Hz. Un latido rapido -- 100 bpm son
% 1.67 Hz -- se situa sobre esa pendiente: aunque forme un pico local
% perfectamente visible, su ALTURA ABSOLUTA es menor que el hombro, y el
% selector, que compara alturas, elige el hombro.
%
% Dividir por una mediana movil convierte cada punto en "cuanto sobresale
% respecto a su entorno". Asi un pico modesto sobre fondo bajo puede ganar
% a un hombro alto sobre fondo alto, que es justo lo que hace falta.
%
% Efecto secundario a vigilar: el blanqueo tambien realza ruido en zonas
% planas. Si empeora casos que antes iban bien, ponlo en false y compara.
blanqueo_hr = false;
if blanqueo_hr
    df  = f_hr_ax(2) - f_hr_ax(1);
    nsm = max(5, 2*round(0.25/df) + 1);      % ventana de ~0.5 Hz, impar
    base = movmedian(P_hr, nsm);
    P_hr = P_hr ./ max(base, max(base)*1e-6);
    fprintf('Blanqueo espectral aplicado (ventana %.2f Hz).\n', nsm*df);
end
% ======================================================
info.f_hr_axis = f_hr_ax;
info.P_hr = P_hr;

delta = max(0.05, 0.10 * f_rr_est);
mask_valid = (f_hr_ax >= fmin_hr) & (f_hr_ax <= fmax_hr);

kill = false(size(f_hr_ax));
if ~apnea
    % Solo tiene sentido enmascarar armonicos si hay respiracion. En apnea,
    % f_rr_est es ruido y enmascarar alrededor de sus multiplos borraria
    % trozos de la banda cardiaca sin motivo.
    for h = 1:3
        kill = kill | (abs(f_hr_ax - h * f_rr_est) <= delta);
    end
end
mask_final = mask_valid & ~kill;
if ~any(mask_final)
    mask_final = mask_valid;
end

if ~any(mask_final)
    info.hr_fallback = true;
    bpm_est = 75;
    info.f_hr_est = 75/60;
    return;
end

[f_hr_est, hr_quality] = select_hr_peak_robust(f_hr_ax, P_hr, mask_final, f_rr_est);
bpm_est = 60 * f_hr_est;
info.f_hr_est = f_hr_est;

info.hr_prominence  = hr_quality.prominence;
info.hr_score_ratio = hr_quality.score_ratio;

% ---- Consistencia por ventanas: la metrica de confianza que de verdad
% ---- vale para una presentacion. Un pico que aparece en una sola FFT de
% ---- 50 s no dice nada; uno que aparece en 8 de 10 ventanas, si.
[info.hr_window_bpm, info.hr_agreement] = ...
    consistencia_hr(x_hr, Fs, hr_win_s, fmin_hr, fmax_hr, bpm_est);

fprintf('\n--- CONSISTENCIA HR ---\n');
fprintf('HR por ventana (bpm): %s\n', formatear_lista(info.hr_window_bpm));
fprintf('Acuerdo con el global (+-5 bpm): %.0f%%\n', 100*info.hr_agreement);
fprintf('-----------------------\n');

info.hr_low_confidence = hr_quality.prominence < 2.0 || ...
                         hr_quality.score_ratio < 1.12 || ...
                         info.hr_agreement < 0.6;

if info.hr_low_confidence
    warning(['HR con baja confianza: prominencia %.2f, ratio %.2f, acuerdo ' ...
        '%.0f%%. Probablemente domina residuo respiratorio, multipath o ' ...
        'movimiento corporal.'], hr_quality.prominence, hr_quality.score_ratio, ...
        100*info.hr_agreement);
end

bpm_est = double(max(30, min(200, bpm_est(1))));
end

% =========================================================================
% AJUSTE DE CIRCULO
% =========================================================================
function [c_t, c0, R0, n_fits] = estimar_centro_circulo(z, Fs, win_s)
%ESTIMAR_CENTRO_CIRCULO Centro del circulo IQ, por ventanas e interpolado.
%   Devuelve c_t (un centro por muestra), mas el ajuste global c0/R0.

z = z(:);
N = numel(z);

[c0, R0] = fit_circle_kasa(z);

Nw = round(win_s * Fs);
hop = max(1, round(Nw/4));
if Nw >= N || Nw < 32
    c_t = repmat(c0, N, 1);
    n_fits = 0;
    return;
end

centros = [];
tiempos = [];
for i0 = 1:hop:(N - Nw + 1)
    seg = z(i0:i0+Nw-1);

    % Solo confiamos en ventanas con arco suficiente: con un arco casi
    % puntual el ajuste de Kasa se dispara hacia el infinito.
    ang = unwrap(angle(seg - c0));
    if (max(ang) - min(ang)) < 0.35
        continue;
    end

    [cw, Rw] = fit_circle_kasa(seg);
    if ~isfinite(cw) || ~isfinite(Rw)
        continue;
    end
    % Guardias de cordura frente al ajuste global.
    if abs(cw - c0) > 0.6 * R0, continue; end
    if Rw < 0.4*R0 || Rw > 2.5*R0, continue; end

    centros(end+1, 1) = cw;      %#ok<AGROW>
    tiempos(end+1, 1) = i0 + (Nw-1)/2;   %#ok<AGROW>
end

n_fits = numel(centros);
if n_fits < 3
    c_t = repmat(c0, N, 1);
    return;
end

% Suavizado ligero y luego interpolacion a cada muestra.
k_smooth = min(5, 2*floor(n_fits/4)+1);
centros = movmean(centros, k_smooth);

idx = (1:N).';
c_t = interp1(tiempos, centros, idx, 'linear');
c_t(idx < tiempos(1))   = centros(1);
c_t(idx > tiempos(end)) = centros(end);
end

function [c, R] = fit_circle_kasa(z)
%FIT_CIRCLE_KASA Ajuste algebraico de circunferencia (Kasa) sobre IQ.
x = real(z(:));
y = imag(z(:));
A = [x, y, ones(numel(x), 1)];
b = x.^2 + y.^2;
p = A \ b;
c = complex(p(1)/2, p(2)/2);
R = sqrt(max(p(3) + abs(c)^2, eps));
end

% =========================================================================
% CANCELACION DE RESPIRACION ADAPTATIVA
% =========================================================================
function x_free = cancelar_respiracion(x, Fs, f_rr0, nharm, win_s)
%CANCELAR_RESPIRACION Ajuste armonico por ventanas deslizantes + overlap-add.
%   En cada ventana se refina f_rr en un +-12% alrededor de f_rr0 y se
%   ajustan nharm armonicos mas una tendencia lineal. El modelo se
%   reconstruye con ventana de Hann y overlap-add, de forma que sigue la
%   deriva natural de la respiracion en vez de imponer una frecuencia fija
%   durante toda la captura.

x = x(:);
N = numel(x);
Nw = round(win_s * Fs);
hop = max(1, round(Nw/2));

if Nw >= N || Nw < 64
    % Captura corta: ajuste unico global.
    x_free = x - modelo_armonico(x, Fs, f_rr0, nharm, 0);
    return;
end

modelo = zeros(N, 1);
pesos  = zeros(N, 1);
w = hann(Nw, 'periodic');

% RANGO DE REFINAMIENTO: +-3%, no +-12%.
%
% Con +-12% un armonico alto se desplaza lo suficiente para aterrizar sobre
% el latido y ser restado como si fuera respiracion. Caso real medido:
% f_rr = 0.2617 Hz, 5o armonico nominal en 1.308 Hz; al refinar a 0.917 se
% mueve a 1.200 Hz, que es exactamente 72 bpm. El modelo lo restaba y el
% pico cardiaco desaparecia.
%
% El margen amplio no hace falta: f_rr se estima con error de decimas de
% rpm (15.7 medido contra 16 real). Con +-3% el 5o armonico solo puede
% moverse entre 1.269 y 1.347 Hz, sin alcanzar la banda del latido.
%
% Regla general: el refinamiento debe cumplir  margen < 1/(2*nharm)  para
% que dos armonicos consecutivos no puedan solaparse. Con nharm = 6 eso da
% un techo del 8%.
refin = 0.03;
f_grid = f_rr0 * linspace(1-refin, 1+refin, 13);

for i0 = 1:hop:(N - Nw + 1)
    seg = x(i0:i0+Nw-1);

    % Refinamiento local de f_rr: nos quedamos con la frecuencia cuyo
    % modelo explica mas energia del segmento.
    mejor_res = Inf;
    mejor_mod = zeros(Nw, 1);
    for ff = f_grid
        m = modelo_armonico(seg, Fs, ff, nharm, i0-1);
        r = sum((seg - m).^2);
        if r < mejor_res
            mejor_res = r;
            mejor_mod = m;
        end
    end

    modelo(i0:i0+Nw-1) = modelo(i0:i0+Nw-1) + w .* mejor_mod;
    pesos(i0:i0+Nw-1)  = pesos(i0:i0+Nw-1) + w;
end

% Los extremos pueden quedar con peso bajo o nulo.
sin_peso = pesos < 1e-6;
modelo(~sin_peso) = modelo(~sin_peso) ./ pesos(~sin_peso);
modelo(sin_peso) = 0;

x_free = x - modelo;
end

function m = modelo_armonico(seg, Fs, f0, nharm, n_offset)
%MODELO_ARMONICO Mejor ajuste por minimos cuadrados de una tendencia lineal
%   mas nharm armonicos de f0 sobre el segmento dado.
Nw = numel(seg);
t = ((0:Nw-1).' + n_offset) / Fs;
t_n = (t - mean(t)) / (max(t) - min(t) + eps);

H = [ones(Nw,1), t_n];
for h = 1:nharm
    H = [H, sin(2*pi*h*f0*t), cos(2*pi*h*f0*t)]; %#ok<AGROW>
end

theta = H \ seg;
m = H * theta;
end

% =========================================================================
% CONSISTENCIA DE HR POR VENTANAS
% =========================================================================
function [bpm_win, acuerdo] = consistencia_hr(x_hr, Fs, win_s, fmin, fmax, bpm_ref)
x_hr = x_hr(:);
N = numel(x_hr);
Nw = round(win_s * Fs);
hop = max(1, round(Nw/2));

bpm_win = [];
if Nw >= N
    bpm_win = bpm_ref;
    acuerdo = 1;
    return;
end

for i0 = 1:hop:(N - Nw + 1)
    seg = x_hr(i0:i0+Nw-1);
    seg = seg - mean(seg);
    Nfft = 2^nextpow2(max(4096, 8*Nw));
    X = abs(fft(seg .* hamming(Nw), Nfft));
    K = floor(Nfft/2) + 1;
    f = (0:K-1).' * Fs / Nfft;
    P = X(1:K);
    m = (f >= fmin) & (f <= fmax);
    if ~any(m), continue; end
    idx = find(m);
    [~, ir] = max(P(idx));
    bpm_win(end+1, 1) = 60 * f(idx(ir)); %#ok<AGROW>
end

if isempty(bpm_win)
    acuerdo = 0;
else
    acuerdo = mean(abs(bpm_win - bpm_ref) <= 5);
end
end

% =========================================================================
% ESPECTRO
% =========================================================================
function [f_ax, P, bin_hz] = espectro(x, Fs, pad_mult, usar_hamming)
if nargin < 4, usar_hamming = true; end
x = x(:);
Nx = numel(x);
Nfft = 2^nextpow2(max(4096, pad_mult * Nx));
bin_hz = Fs / Nfft;
if usar_hamming, wv = hamming(Nx); else, wv = ones(Nx, 1); end
X = fft((x - mean(x)) .* wv, Nfft);
K = floor(Nfft/2) + 1;
f_ax = (0:K-1).' * bin_hz;
P = abs(X(1:K));
end

function txt = formatear_lista(v)
if isempty(v)
    txt = 'sin datos';
    return;
end
parts = strings(1, numel(v));
for i = 1:numel(v)
    parts(i) = sprintf('%.1f', v(i));
end
txt = strjoin(parts, ', ');
end

% =========================================================================
% SELECTORES ESPECTRALES (sin cambios de fondo respecto a v3)
% =========================================================================
function f_est = select_rr_peak_robust(f_ax, P, mask, bin_hz)
idx_band = find(mask);
if isempty(idx_band)
    f_est = NaN;
    return;
end

idx_peaks = local_peak_indices(P, idx_band);
if isempty(idx_peaks)
    [~, i_rel] = max(P(idx_band));
    idx_peaks = idx_band(i_rel);
end

scores = zeros(size(idx_peaks));
f_cands = zeros(size(idx_peaks));
for i = 1:numel(idx_peaks)
    k = idx_peaks(i);
    f0 = parabolic_peak_frequency(f_ax, P, k, bin_hz);
    local_floor = local_noise_floor(f_ax, P, f0, 0.12, 0.025, mask);
    prom = P(k) / (local_floor + eps);

    edge_penalty = 1.0;
    if f0 < 8/60
        edge_penalty = 0.25;
    elseif f0 < 10/60
        edge_penalty = 0.55;
    elseif f0 < 12/60
        edge_penalty = 0.85;
    end

    scores(i) = prom * edge_penalty;
    f_cands(i) = f0;
end

scores(~isfinite(scores)) = 0;
if ~any(scores > 0)
    [~, i_rel] = max(P(idx_band));
    f_est = f_ax(idx_band(i_rel));
else
    [~, best] = max(scores);
    f_est = f_cands(best);
end

f_est = prefer_rr_fundamental(f_est, f_ax, P, mask, bin_hz);
end

function f_est = prefer_rr_fundamental(f_est, f_ax, P, mask, bin_hz)
if ~isfinite(f_est) || f_est < 24/60
    return;
end

k_est = nearest_masked_bin(f_ax, mask, f_est);
if isempty(k_est)
    return;
end
p_est = P(k_est);

best_f = f_est;
best_ratio = 0;
for divisor = [2 3]
    f_candidate = f_est / divisor;
    if f_candidate < 8/60 || f_candidate > 24/60
        continue;
    end
    if divisor == 3 && (f_est < 33/60 || f_candidate < 10/60)
        continue;
    end

    k0 = nearest_masked_bin(f_ax, mask, f_candidate);
    if isempty(k0)
        continue;
    end

    half_width_bins = max(2, round(0.035 / bin_hz));
    lo = max(1, k0 - half_width_bins);
    hi = min(numel(P), k0 + half_width_bins);
    local_idx = lo:hi;
    local_idx = local_idx(mask(local_idx));
    if isempty(local_idx)
        continue;
    end

    [p_fund, i_rel] = max(P(local_idx));
    f_fund = f_ax(local_idx(i_rel));
    ratio = p_fund / (p_est + eps);
    required_ratio = 0.45;
    if divisor == 3
        required_ratio = 0.55;
    end

    if ratio >= required_ratio && ratio > best_ratio
        best_ratio = ratio;
        best_f = f_fund;
    end
end

if best_ratio > 0
    f_est = best_f;
end
end

function k = nearest_masked_bin(f_ax, mask, f0)
idx = find(mask);
if isempty(idx)
    k = [];
    return;
end
[~, i_rel] = min(abs(f_ax(idx) - f0));
k = idx(i_rel);
end

function [f_est, quality] = select_hr_peak_robust(f_ax, P, mask, f_rr)
quality = struct('prominence', 0, 'score', 0, 'score_ratio', 0);
idx_band = find(mask);
if isempty(idx_band)
    f_est = 75 / 60;
    return;
end

idx_peaks = local_peak_indices(P, idx_band);
if isempty(idx_peaks)
    [~, i_rel] = max(P(idx_band));
    idx_peaks = idx_band(i_rel);
end

scores = zeros(size(idx_peaks));
proms  = zeros(size(idx_peaks));
f_cands = zeros(size(idx_peaks));
for i = 1:numel(idx_peaks)
    k = idx_peaks(i);
    f0 = f_ax(k);
    local_floor = local_noise_floor(f_ax, P, f0, 0.35, 0.06, mask);
    prom = P(k) / (local_floor + eps);
    proms(i) = prom;

    harm_pen = 0;
    delta_low  = max(0.05,  0.10 * f_rr);
    delta_high = max(0.035, 0.07 * f_rr);
    for h = 1:3
        harm_pen = harm_pen + 1.0 * exp(-((f0 - h*f_rr).^2) / (2 * delta_low^2));
    end
    for h = 4:6
        harm_pen = harm_pen + 0.12 * exp(-((f0 - h*f_rr).^2) / (2 * delta_high^2));
    end

    adult_bonus = 1.0;
    if f0 >= 60/60 && f0 <= 130/60
        adult_bonus = 1.15;
    end

    high_rate_penalty = 1.0;
    if f0 > 130/60
        high_rate_penalty = 0.70;
    end
    if f0 > 150/60
        high_rate_penalty = 0.55;
    end

    scores(i) = prom * adult_bonus * high_rate_penalty / (1 + 1 * harm_pen);
    f_cands(i) = f0;
end

scores(~isfinite(scores)) = 0;
if ~any(scores > 0)
    [~, i_rel] = max(P(idx_band));
    f_est = f_ax(idx_band(i_rel));
else
    [~, best] = max(scores);
    f_est = f_cands(best);
    quality.prominence = proms(best);
    quality.score      = scores(best);
    sorted_scores = sort(scores(:), 'descend');
    if numel(sorted_scores) >= 2 && sorted_scores(2) > 0
        quality.score_ratio = sorted_scores(1) / sorted_scores(2);
    else
        quality.score_ratio = Inf;
    end
end
end

function idx_peaks = local_peak_indices(P, idx_band)
idx_peaks = [];
for ii = 1:numel(idx_band)
    k = idx_band(ii);
    if k > 1 && k < numel(P) && P(k) >= P(k-1) && P(k) >= P(k+1)
        idx_peaks(end+1) = k; 
    end
end
end

function f0 = parabolic_peak_frequency(f_ax, P, k, bin_hz)
f0 = f_ax(k);
if k <= 1 || k >= numel(P)
    return;
end
al = P(k-1); be = P(k); ga = P(k+1);
denom = al - 2*be + ga;
if abs(denom) > 1e-12
    delta_idx = 0.5 * (al - ga) / denom;
    delta_idx = max(-0.5, min(0.5, delta_idx));
    f0 = f_ax(k) + delta_idx * bin_hz;
end
end

function nf = local_noise_floor(f_ax, P, f0, half_width, guard_width, mask)
idx = mask & abs(f_ax - f0) <= half_width & abs(f_ax - f0) >= guard_width;
if nnz(idx) < 8
    idx = mask;
end
nf = median(P(idx));
if ~isfinite(nf) || nf <= 0
    nf = median(P(mask)) + eps;
end
end

function check_filter_fs(d, name, fs_runtime)
try
    fs_design = d.SampleRate;
catch
    return;
end
if abs(fs_design - fs_runtime) > 1e-6
    warning(['%s fue disenado para Fs=%.3f Hz pero se esta llamando con ' ...
        'Fs=%.3f Hz. El filtro actuara sobre la banda equivocada sin ' ...
        'lanzar error.'], name, fs_design, fs_runtime);
end
end

function [z_out, f_off, vueltas, deriva_pp] = quitar_offset_portadora(z, Fs, fc_track)
%QUITAR_OFFSET_PORTADORA Elimina la rotacion residual de fase por offset de LO.
%
%   Distingue el offset de LO de la respiracion por una propiedad
%   estructural, no por frecuencia: la respiracion es una oscilacion
%   ACOTADA (la fase vuelve siempre a su linea base), mientras que un
%   offset de portadora ACUMULA fase monotonamente. Un ajuste de primer
%   orden sobre la fase desenrollada da pendiente ~0 para una oscilacion
%   acotada, por muy lenta que sea, y la pendiente real para una rampa.
%
%   Por eso es seguro aunque el offset (aqui ~0.1 Hz) caiga dentro de la
%   banda respiratoria: no se filtra por frecuencia, se separa por forma.

if nargin < 3 || isempty(fc_track), fc_track = 0.06; end   % Hz (0 = desactivado)

N = numel(z);
n = (0:N-1).';

% --- Etapa 1: pico espectral grueso ------------------------------------
% La componente mas fuerte del baseband es el clutter estatico + la fuga
% TX->RX. Sin offset deberia estar exactamente en DC; donde este es el
% offset.
Nf = 2^nextpow2(8*N);
Zf = fft(z, Nf);
f_ejes = (0:Nf-1).' / Nf * Fs;
f_ejes(f_ejes > Fs/2) = f_ejes(f_ejes > Fs/2) - Fs;

mag = abs(Zf);
mag(abs(f_ejes) > 2) = 0;        % un offset de LO no pasa de unos pocos Hz
[~, k] = max(mag);
f0 = f_ejes(k);

z1 = z .* exp(-1j*2*pi*f0*n/Fs);

% --- Etapa 2: refinamiento por ajuste lineal de fase --------------------
psi = unwrap(angle(z1));
pcoef = polyfit(n, psi, 1);
f1 = pcoef(1) * Fs / (2*pi);

z2 = z1 .* exp(-1j*2*pi*f1*n/Fs);

f_off   = f0 + f1;
vueltas = f_off * N / Fs;

% --- Etapa 3: seguimiento de la DERIVA lenta ---------------------------
% Las etapas 1 y 2 quitan un offset CONSTANTE. Pero los PLL de TX y RX son
% independientes y su ruido de fase cercano a portadora no esta
% correlacionado: integrado sobre 60 s, la diferencia de fase hace un
% paseo lento de varios radianes que ninguna recta puede seguir.
%
% Se estima esa deriva como la parte de la fase por debajo de fc_track y se
% derrota con ella. La separacion sigue siendo estructural: la respiracion
% de un adulto en reposo esta en 0.20-0.35 Hz (12-21 rpm), muy por encima
% del corte. Con respiracion MUY lenta (< 8 rpm) este paso si mordería
% senal util; por eso se imprime la deriva retirada, para poder vigilarlo.
if fc_track <= 0
    z_out     = z2;          % desactivado (captura corta)
    deriva_pp = 0;
else
    psi = unwrap(angle(z2));
    Wn  = fc_track / (Fs/2);
    [bb, aa] = butter(2, Wn, 'low');
    psi_lento = filtfilt(bb, aa, psi);

    z_out    = z2 .* exp(-1j*psi_lento);
    deriva_pp = max(psi_lento) - min(psi_lento);
end
end


function [f_ax, P] = espectro_welch(x, Fs, win_s, hop_frac, pad_mult)
%ESPECTRO_WELCH Periodograma promediado sobre ventanas solapadas.
%
%   win_s     longitud de ventana en segundos
%   hop_frac  avance como fraccion de la ventana (0.25 = 75% de solape)
%
%   Si la senal no da ni para una ventana, cae de vuelta a una FFT unica.

N = numel(x);
L = round(win_s * Fs);

if L >= N
    [f_ax, P, ~] = espectro(x, Fs, pad_mult);
    return;
end

hop = max(1, round(hop_frac * L));
w   = hann(L);
w   = w / sqrt(mean(w.^2));        % conserva la potencia

inicios = 1:hop:(N - L + 1);

Pmat = [];
for i = 1:numel(inicios)
    seg = x(inicios(i) : inicios(i)+L-1);
    seg = seg - mean(seg);
    [f_ax, Pi, ~] = espectro(seg .* w, Fs, pad_mult, false);
    Pmat(:, end+1) = Pi(:);
end

P = mean(Pmat, 2);
end
