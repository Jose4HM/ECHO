# Guía de trazabilidad: adquisición I/Q y diagnóstico RF de ECHO

**Archivo de referencia:** `vacio01.mat`\
**Objetivo:** documentar de dónde procede cada valor utilizado en las
subsecciones *Receiver amplitude and saturation assessment* y *RF
acquisition characterization*, cómo se calcula y cómo reproducirlo en
MATLAB.

> **Criterio de trazabilidad:** se distingue entre (1) valores
> calculados directamente a partir de `rx_continuous`, (2) valores
> leídos de las estructuras guardadas en el `.mat` y (3)
> interpretaciones que requieren comprobar el código de adquisición. No
> se deben presentar como mediciones físicas calibradas las amplitudes
> digitales I/Q.

------------------------------------------------------------------------

## 1. Archivo y variables relevantes

El archivo `vacio01.mat` contiene, entre otras, estas variables:

  ---------------------------------------------------------------------------
  Variable                    Dimensiones Tipo               Uso
  ----------------- --------------------- ------------------ ----------------
  `rx_continuous`          60,000,000 × 1 `double complex`   Muestras I/Q
                                                             continuas
                                                             analizadas

  `Fs_frontend`                     1 × 1 `double`           Frecuencia de
                                                             muestreo
                                                             registrada

  `cfg`                             1 × 1 `struct`           Configuración y
                                                             diagnóstico de
                                                             adquisición

  `vital_info`                      1 × 1 `struct`           Diagnósticos del
                                                             procesamiento
                                                             fisiológico

  `bpm`                             1 × 1 `double`           Salida HR
                                                             guardada en el
                                                             archivo

  `rpm`                             1 × 1 `double`           Salida RR
                                                             guardada en el
                                                             archivo

  `iq_frontend`                 6,000 × 1 `double complex`   Otra señal
                                                             almacenada; no
                                                             es la señal
                                                             continua de 60
                                                             millones de
                                                             muestras
  ---------------------------------------------------------------------------

La variable principal para reproducir las estadísticas de amplitud es
`rx_continuous`. Tiene 60 millones de muestras complejas y ocupa
aproximadamente 916 MB. Por eso, el script final la procesa por bloques
mediante `matfile`, sin cargarla completa en memoria.

## 2. Duración y cantidad de muestras

La configuración guardada indica:

-   Frecuencia de muestreo: (F_s=1,000,000 `\text{Hz}`{=tex}) (1 MS/s).
-   Duración configurada: (T=60 `\text{s}`{=tex}).
-   Número de bloques: 6,000.
-   Longitud de bloque: 10,000 muestras.
-   Muestras totales registradas: 60,000,000.

La cantidad esperada se calcula como:

\[ N\_{`\mathrm{esperado}`{=tex}}=F_sT \]

\[ N\_{`\mathrm{esperado}`{=tex}}=1,000,000`\times60`{=tex}
=60,000,000 `\text{muestras}`{=tex} \]

También coincide con el producto del número de bloques por la longitud
de cada bloque:

\[ N=6000`\times10000`{=tex}=60,000,000 \]

En `vacio01.mat`, `cfg.total_actual_samples` registra 60,000,000 y
`rx_continuous` tiene dimensiones 60,000,000 × 1. La coincidencia
confirma que la cantidad almacenada concuerda con la duración y
frecuencia de muestreo configuradas. Por sí sola, la cantidad de
muestras no demuestra ausencia de pérdidas; para ello se consultan
también los contadores de integridad guardados.

## 3. Representación de las muestras I/Q

Cada muestra compleja se representa como:

\[ x\[n\]=I\[n\]+jQ\[n\] \]

donde (I\[n\]) es la componente en fase, (Q\[n\]) es la componente en
cuadratura, (j) es la unidad imaginaria y (n) es el índice de muestra.

La magnitud de una muestra compleja es:

\[ \|x\[n\]\|=`\sqrt{I[n]^2+Q[n]^2}`{=tex} \]

Esta magnitud es una amplitud digital. No es, por sí sola, una potencia
RF en vatios ni una potencia recibida en dBm.

## 4. Máxima magnitud compleja

La máxima magnitud se obtiene buscando el mayor valor de (\|x\[n\]\|) en
las 60 millones de muestras:

\[
A\_{`\max`{=tex}}=`\max`{=tex}\_{1`\leq `{=tex}n`\leq `{=tex}N}\|x\[n\]\|
\]

En MATLAB, la operación equivalente es:

``` matlab
max(abs(rx_continuous))
```

El procesamiento por bloques reprodujo:

\[ `\boxed{A_{\max}=0.947688875}`{=tex} \]

Redondeado a seis decimales:

\[ A\_{`\max`{=tex}}=0.947689 \]

Este es el valor que aparece en el paper. Además, `cfg.peak_global`
guarda aproximadamente 0.9477; el cálculo independiente desde
`rx_continuous` permite verificarlo con mayor precisión.

## 5. Amplitud RMS compleja

La amplitud RMS compleja se calcula como la raíz cuadrada del promedio
de los cuadrados de las magnitudes complejas:

\[ A\_{`\mathrm{RMS}`{=tex}}=
`\sqrt{\frac{1}{N}\sum_{n=1}^{N}|x[n]|^2}`{=tex} \]

Como (\|x\[n\]\|^2=I\[n\]^2+Q\[n\]\^2), también puede escribirse:

\[ A\_{`\mathrm{RMS}`{=tex}}= `\sqrt{\frac{1}{N}\sum_{n=1}^{N}
\left(I[n]^2+Q[n]^2\right)}`{=tex} \]

En MATLAB:

``` matlab
sqrt(mean(abs(rx_continuous).^2))
```

El cálculo independiente por bloques produjo:

\[ `\boxed{A_{\mathrm{RMS}}=0.465224451}`{=tex} \]

Redondeado a seis decimales:

\[ A\_{`\mathrm{RMS}`{=tex}}=0.465224 \]

Este valor describe el nivel cuadrático medio de las muestras digitales
almacenadas. No debe denominarse potencia RF calibrada ni convertirse
directamente a vatios sin una calibración de la cadena receptora.

## 6. Comprobación de saturación I/Q

En esta evaluación se utilizó un umbral de saturación definido por
software:

\[ A\_{`\mathrm{umbral}`{=tex}}=0.99 \]

Se revisan por separado los valores absolutos de las componentes I y Q:

\[ I\_{`\max`{=tex}}=`\max`{=tex}\_n\|I\[n\]\| \]

\[ Q\_{`\max`{=tex}}=`\max`{=tex}\_n\|Q\[n\]\| \]

La condición de ausencia de muestras que alcancen el umbral es:

\[ I\_{`\max`{=tex}}\<0.99,`\qquad `{=tex}Q\_{`\max`{=tex}}\<0.99 \]

Los resultados reproducidos desde `rx_continuous` fueron:

  Magnitud                   Resultado              Umbral
  ---------------------- ------------- -------------------
  Máximo absoluto de I     0.521484375                0.99
  Máximo absoluto de Q     0.947265625                0.99
  Muestras I con (                   I   `\geq0.99`{=tex})
  Muestras Q con (                   Q   `\geq0.99`{=tex})

Por tanto, ninguna muestra I ni Q alcanzó el umbral definido por
software.

El máximo complejo también fue inferior a 0.99 (0.947688875). Como
(\|I\[n\]\|`\leq`{=tex}\|x\[n\]\|) y (\|Q\[n\]\|`\leq`{=tex}\|x\[n\]\|),
esto es coherente con los máximos individuales. De todos modos, el
script cuenta explícitamente las muestras I y Q que alcanzan el umbral.

**Limitación de interpretación:** 0.99 es el umbral aplicado por el
software para esta comprobación. No debe describirse automáticamente
como el límite físico exacto del ADC ni como una prueba de ausencia de
cualquier tipo de distorsión.

## 7. Integridad de la adquisición

Los siguientes valores no se calculan a partir de la magnitud de
`rx_continuous`; se leen de la estructura `cfg` guardada durante la
adquisición:

  Campo de `cfg`             Valor guardado Interpretación
  ------------------------ ---------------- -------------------------------------------------
  `overruns`                              0 No se registraron buffer overruns
  `n_incomplete_blocks`                   0 No se registraron bloques incompletos
  `total_actual_samples`         60,000,000 Total de muestras registrado
  `capture_valid`                         1 La captura fue marcada como válida
  `saturation`                            0 Indicador de saturación guardado por el proceso

Estos campos sustentan la frase del paper sobre la integridad de la
adquisición. Es importante atribuirlos al registro de estado de la
captura, no a un cálculo posterior de las muestras.

## 8. Comparación de las condiciones baseline

La estructura `cfg` también guarda las amplitudes RMS de las condiciones
con el transmisor deshabilitado y habilitado:

-   `cfg.baseline_rms_off`: RMS con TX deshabilitado.
-   `cfg.baseline_rms_on`: RMS con TX habilitado.

Los valores mostrados en MATLAB fueron aproximadamente:

\[ A\_{`\mathrm{off}`{=tex}}=0.0091 \]

\[ A\_{`\mathrm{on}`{=tex}}=0.4764 \]

La diferencia relativa de amplitud en decibelios se calcula mediante:

\[ `\Delta`{=tex}*{`\mathrm{dB}`{=tex}}=
20`\log`{=tex}*{10}`\left`{=tex}(
`\frac{A_{\mathrm{on}}}{A_{\mathrm{off}}}`{=tex} `\right`{=tex}) \]

Se utiliza (20`\log`{=tex}\_{10}) porque se comparan amplitudes RMS, no
potencias.

El campo `cfg.leak_over_floor_db` guarda el valor reportado por el
procedimiento de adquisición: 34.3690 dB. **No se debe intentar
reproducir exactamente ese valor usando únicamente los valores 0.0091 y
0.4764 redondeados a cuatro decimales.** Al hacerlo, se obtiene
aproximadamente 34.3786 dB. La diferencia puede deberse a que los campos
RMS internos conservan más decimales que los mostrados por `disp`, o a
la precisión de los valores usados por el cálculo original. El script
final imprime los valores con mayor precisión y compara el cálculo con
el campo guardado.

La diferencia entre las condiciones TX apagado y encendido incluye el
acoplamiento directo Tx--Rx y las reflexiones del entorno. No representa
por sí sola el aislamiento entre antenas.

## 9. Diagnósticos del procesamiento fisiológico

Los diagnósticos se leen de `vital_info`. No se recalculan a partir de
`rx_continuous` en este procedimiento.

  --------------------------------------------------------------------------
  Campo                                Valor reportado Significado
  ----------------------- ---------------------------- ---------------------
  `carrier_offset_hz`                       −0.0554 Hz Desplazamiento de
                                                       frecuencia de
                                                       portadora estimado

  `arc_span_rad`                            3.8911 rad Extensión angular de
                                                       la trayectoria I/Q

  `center_offset_ratio`                         0.7172 Relación normalizada
                                                       de desplazamiento del
                                                       centro

  `circle_radius`                               0.2751 Radio del ajuste
                                                       circular en la escala
                                                       de la señal

  `rr_candidate_rpm`                       19.1534 RPM Candidato de
                                                       frecuencia
                                                       respiratoria

  `hr_candidate_bpm`                      129.2725 BPM Candidato de
                                                       frecuencia cardíaca

  `rr_low_confidence`                                1 Indicador de baja
                                                       confianza para RR

  `hr_low_confidence`                                1 Indicador de baja
                                                       confianza para HR

  `rr_window_agreement`                         0.5000 Acuerdo entre
                                                       ventanas para RR

  `hr_agreement`                                0.2500 Acuerdo entre
                                                       ventanas para HR
  --------------------------------------------------------------------------

Los campos `rr_candidate_rpm` y `hr_candidate_bpm` son candidatos, no
necesariamente estimaciones aceptadas. Para afirmar que las salidas
finales fueron `NaN`, hay que comprobar los campos superiores `rpm` y
`bpm` almacenados en el archivo. El script final los muestra por
separado junto con los candidatos y los indicadores de baja confianza.

En esta captura, ambos indicadores de baja confianza valen 1. La
aceptación o el rechazo final debe describirse según las salidas
efectivamente guardadas y las reglas de decisión del algoritmo.

## 10. Limitaciones de la configuración de ganancia RX

La configuración guardada incluye:

-   `cfg.rx_channel = 'RX2'`.
-   `cfg.rx_gain = 0`.
-   `cfg.rx_agc = 'manual'`.
-   `cfg.rx_gain_actual = 0`.
-   `cfg.rx_gain_mode_prop = 'agc'`.

La inspección del archivo `bladeRF_XCVR.m` mostró que los métodos
`set.agc`, `get.agc`, `set.gain` y `get.gain` seleccionan explícitamente
`BLADERF_CHANNEL_RX1` cuando el objeto es de dirección RX. En cambio, la
adquisición de ECHO se realizó por RX2.

Por ello, las llamadas de configuración y consulta de ganancia de ese
binding no controlan ni verifican RX2 mediante esos métodos. El valor
`rx_gain=0` no debe presentarse como la ganancia efectiva aplicada a
RX2. Asimismo, `rx_gain_mode_prop='agc'` procede de una consulta cuyo
método está fijado a RX1; no demuestra por sí solo el modo efectivo de
RX2.

Se ha indicado que el receptor termina usando su configuración
predeterminada, descrita como 60 dB, pero ese valor exacto no queda
demostrado por las estadísticas de `vacio01.mat` ni por los campos
anteriores. Para afirmarlo como ganancia efectiva se necesita una
lectura independiente del canal RX2 o documentación técnica que
establezca inequívocamente ese estado.

Esta limitación no invalida las estadísticas calculadas sobre las
muestras digitales guardadas. Sí debe considerarse al interpretar
amplitudes y comparar adquisiciones.

------------------------------------------------------------------------

## 11. Script MATLAB completo de reproducción

**Instrucciones:** guarda este código como `verificar_vacio01_IQ.m` en
la misma carpeta que `vacio01.mat` y ejecútalo desde MATLAB. El script
utiliza `matfile` y procesa la señal por bloques de un millón de
muestras, evitando cargar los aproximadamente 916 MB de `rx_continuous`
de una sola vez.

El código separa los resultados calculados directamente desde las
muestras de los campos leídos de `cfg` y `vital_info`.

``` matlab
%% VERIFICACION DE TRAZABILIDAD I/Q - ECHO
% Archivo: vacio01.mat
% Procesa rx_continuous por bloques para evitar cargar
% toda la señal compleja en memoria.

clear;
clc;

archivo = 'vacio01.mat';

if ~isfile(archivo)
    error('No se encuentra el archivo: %s', archivo);
end

%% 1. Cargar metadatos y estructuras guardadas

S = load(archivo, ...
    'Fs_frontend', 'cfg', 'vital_info', 'bpm', 'rpm');

cfg = S.cfg;
vital_info = S.vital_info;

M = matfile(archivo);

% Obtener dimensiones sin cargar rx_continuous completo
info = whos('-file', archivo);
idx = find(strcmp({info.name}, 'rx_continuous'), 1);

if isempty(idx)
    error('No existe rx_continuous en el archivo.');
end

dims = info(idx).size;

if numel(dims) ~= 2 || dims(2) ~= 1
    error('Se esperaba rx_continuous como vector columna.');
end

N_real = dims(1);

%% 2. Configuracion temporal

Fs = cfg.fs;
T = cfg.capture_seconds;
N_esperado = Fs * T;

bloque = 1e6;
umbral = 0.99;

fprintf('\n========================================\n');
fprintf('1. CONFIGURACION Y CANTIDAD DE MUESTRAS\n');
fprintf('========================================\n');

fprintf('Fs en cfg:              %.0f Hz\n', Fs);
fprintf('Fs_frontend:            %.0f Hz\n', S.Fs_frontend);
fprintf('Duracion configurada:   %.6f s\n', T);
fprintf('Muestras esperadas:     %.0f\n', N_esperado);
fprintf('Muestras en archivo:    %.0f\n', N_real);
fprintf('Muestras en cfg:        %.0f\n', ...
    cfg.total_actual_samples);
fprintf('Bloques configurados:   %d\n', cfg.num_blocks);
fprintf('Longitud de bloque:     %d\n', cfg.block_len);

%% 3. Inicializar acumuladores

suma_cuadrados = 0;

max_complejo = 0;
max_I = 0;
max_Q = 0;

cuenta_sat_I = 0;
cuenta_sat_Q = 0;

%% 4. Procesar rx_continuous por bloques

fprintf('\n========================================\n');
fprintf('2. PROCESAMIENTO DE MUESTRAS I/Q\n');
fprintf('========================================\n');

for inicio = 1:bloque:N_real

    fin = min(inicio + bloque - 1, N_real);

    % MatFile requiere indicar las dos dimensiones.
    x = M.rx_continuous(inicio:fin, 1);

    I = real(x);
    Q = imag(x);

    magnitud = abs(x);

    % RMS: acumulacion de cuadrados de magnitud
    suma_cuadrados = suma_cuadrados + ...
        sum(magnitud.^2);

    % Maximos globales
    max_complejo = max(max_complejo, max(magnitud));
    max_I = max(max_I, max(abs(I)));
    max_Q = max(max_Q, max(abs(Q)));

    % Conteo de muestras que alcanzan el umbral
    cuenta_sat_I = cuenta_sat_I + ...
        nnz(abs(I) >= umbral);

    cuenta_sat_Q = cuenta_sat_Q + ...
        nnz(abs(Q) >= umbral);

    fprintf('Procesadas: %d / %d\n', fin, N_real);
end

% RMS global
rms_complejo = sqrt(suma_cuadrados / N_real);

fprintf('\n========================================\n');
fprintf('3. AMPLITUD Y SATURACION CALCULADAS\n');
fprintf('========================================\n');

fprintf('Maximo complejo: %.12f\n', max_complejo);
fprintf('RMS complejo:    %.12f\n', rms_complejo);
fprintf('Maximo |I|:      %.12f\n', max_I);
fprintf('Maximo |Q|:      %.12f\n', max_Q);

fprintf('\nUmbral software: %.6f\n', umbral);
fprintf('Muestras I >= umbral: %d\n', cuenta_sat_I);
fprintf('Muestras Q >= umbral: %d\n', cuenta_sat_Q);

%% 5. Leer indicadores de integridad guardados

fprintf('\n========================================\n');
fprintf('4. INTEGRIDAD GUARDADA EN CFG\n');
fprintf('========================================\n');

fprintf('Overruns:             %d\n', cfg.overruns);
fprintf('Bloques incompletos:  %d\n', ...
    cfg.n_incomplete_blocks);
fprintf('Total actual samples: %.0f\n', ...
    cfg.total_actual_samples);
fprintf('Capture valid:        %d\n', cfg.capture_valid);
fprintf('Saturation flag:      %d\n', cfg.saturation);

%% 6. Baseline con TX apagado y encendido

fprintf('\n========================================\n');
fprintf('5. BASELINE TX OFF / TX ON\n');
fprintf('========================================\n');

rms_off = cfg.baseline_rms_off;
rms_on = cfg.baseline_rms_on;

delta_baseline_db = 20 * log10(rms_on / rms_off);

fprintf('Baseline RMS OFF: %.12f\n', rms_off);
fprintf('Baseline RMS ON:  %.12f\n', rms_on);
fprintf('Incremento calculado: %.8f dB\n', ...
    delta_baseline_db);

fprintf('Incremento guardado cfg.leak_over_floor_db: %.8f dB\n', ...
    cfg.leak_over_floor_db);

%% 7. Diagnosticos fisiologicos guardados

fprintf('\n========================================\n');
fprintf('6. DIAGNOSTICOS VITAL_INFO\n');
fprintf('========================================\n');

fprintf('Carrier offset:       %.9f Hz\n', ...
    vital_info.carrier_offset_hz);
fprintf('Arc span:             %.9f rad\n', ...
    vital_info.arc_span_rad);
fprintf('Center offset ratio:  %.9f\n', ...
    vital_info.center_offset_ratio);
fprintf('Circle radius:        %.9f\n', ...
    vital_info.circle_radius);

fprintf('RR candidate:         %.9f RPM\n', ...
    vital_info.rr_candidate_rpm);
fprintf('HR candidate:         %.9f BPM\n', ...
    vital_info.hr_candidate_bpm);

fprintf('RR low confidence:    %d\n', ...
    vital_info.rr_low_confidence);
fprintf('HR low confidence:    %d\n', ...
    vital_info.hr_low_confidence);

fprintf('RR window agreement:  %.9f\n', ...
    vital_info.rr_window_agreement);
fprintf('HR agreement:         %.9f\n', ...
    vital_info.hr_agreement);

%% 8. Salidas finales guardadas en el archivo

fprintf('\n========================================\n');
fprintf('7. SALIDAS FINALES GUARDADAS\n');
fprintf('========================================\n');

fprintf('bpm (salida HR): %.9f\n', S.bpm);
fprintf('rpm (salida RR): %.9f\n', S.rpm);

%% 9. Configuracion de ganancia y canales

fprintf('\n========================================\n');
fprintf('8. CONFIGURACION RX / TX\n');
fprintf('========================================\n');

fprintf('TX channel:          %s\n', cfg.tx_channel);
fprintf('RX channel:          %s\n', cfg.rx_channel);
fprintf('TX gain solicitado:  %.6f dB\n', cfg.tx_gain);
fprintf('RX gain solicitado:  %.6f dB\n', cfg.rx_gain);
fprintf('RX gain actual:      %.6f dB\n', cfg.rx_gain_actual);
fprintf('RX AGC solicitado:   %s\n', cfg.rx_agc);
fprintf('RX mode property:    %s\n', cfg.rx_gain_mode_prop);

%% 10. Comparacion con los valores del paper

fprintf('\n========================================\n');
fprintf('9. COMPARACION CON EL PAPER\n');
fprintf('========================================\n');

fprintf('Maximo paper: 0.947689\n');
fprintf('Maximo MATLAB: %.12f\n', max_complejo);

fprintf('RMS paper: 0.465224\n');
fprintf('RMS MATLAB: %.12f\n', rms_complejo);

fprintf('\n===== FIN DE LA VERIFICACION =====\n');
```

------------------------------------------------------------------------

## 12. Resultados que ya se han reproducido

La ejecución realizada sobre `vacio01.mat` produjo:

  Parámetro                           Resultado reproducido
  --------------------------------- -----------------------
  Máxima magnitud compleja                      0.947688875
  RMS complejo                                  0.465224451
  Máximo absoluto I                             0.521484375
  Máximo absoluto Q                             0.947265625
  Muestras I que alcanzan 0.99                            0
  Muestras Q que alcanzan 0.99                            0
  Buffer overruns registrados                             0
  Bloques incompletos registrados                         0
  Total de muestras                              60,000,000
  Captura válida                                          1

Estos resultados son reproducibles a partir de `rx_continuous` y los
registros de `cfg` de `vacio01.mat`.

**Pendientes de verificación adicional:** el modo y la ganancia
efectivos de RX2, el valor exacto de la ganancia predeterminada y la
correspondencia precisa entre `cfg.leak_over_floor_db` y los valores RMS
internos sin redondear.


