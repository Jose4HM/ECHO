% =========================================================================
% ANIMACIÓN DEL CÁLCULO DEL MAE (Mean Absolute Error)
% Datos: Primeras 10 muestras de HR de Joss
% =========================================================================

% 1. Definición de los datos
y_real = [99, 100, 84, 89, 80, 85.5, 85, 88, 88, 97];
y_script = [94.5, 98.5, 95.9, 79.1, 73.2, 82.4, 83.9, 85.7, 89, 98.1];
N = length(y_real);

% 2. Configuración inicial de la ventana (Figura)
figure('Name', 'Animacion MAE', 'Color', 'w', 'Position', [100, 100, 850, 550]);
hold on; grid on;
xlim([0 N+1]);
ylim([min([y_real, y_script]) - 5, max([y_real, y_script]) + 5]);
xlabel('Número de Muestra', 'FontWeight', 'bold');
ylabel('Latidos por minuto (HR)', 'FontWeight', 'bold');

% 3. Dibujar las curvas completas primero
plot(1:N, y_real, '-ob', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'DisplayName', 'HR Real (Azul)');
plot(1:N, y_script, '-sr', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'HR Script (Rojo)');
legend('Location', 'southwest', 'FontSize', 11);

% Título inicial
title('Preparando cálculo del MAE...', 'FontSize', 14);
pause(2); % Pausa de 2 segundos antes de empezar la animación

% 4. Variables para acumular el error
suma_errores = 0;

% 5. Bucle de animación (Paso a paso)
for i = 1:N
    % Calcular la distancia (Error Absoluto) entre el punto real y el script
    error_abs = abs(y_real(i) - y_script(i));

    % Sumar al acumulador total
    suma_errores = suma_errores + error_abs;

    % Calcular el MAE hasta el momento (Promedio)
    mae_parcial = suma_errores / i;

    % --- Animación en la gráfica ---

    % Dibujar una barra negra gruesa que conecte ambos puntos (El error físico)
    plot([i, i], [y_real(i), y_script(i)], 'k', 'LineWidth', 2.5, 'HandleVisibility', 'off');

    % Resaltar los puntos que se están evaluando
    plot(i, y_real(i), 'ob', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(i, y_script(i), 'sr', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');

    % Escribir el valor del error de esa muestra al lado de la línea
    punto_medio_y = (y_real(i) + y_script(i)) / 2;
    text(i + 0.15, punto_medio_y, sprintf('e = %.1f', error_abs), ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');

    % Actualizar el título mostrando las matemáticas en vivo
    title(sprintf('Paso %d/%d | Error de hoy: %.1f | Suma Total: %.1f | MAE: %.2f', ...
        i, N, error_abs, suma_errores, mae_parcial), 'FontSize', 14, 'Color', 'k');

    % Pausa de 1.5 segundos para que tus ojos vean el progreso
    pause(1.5);
end

% 6. Conclusión
title(sprintf('¡Cálculo Finalizado! MAE Final = %.2f', suma_errores / N), ...
    'FontSize', 16, 'Color', [0, 0.5, 0], 'FontWeight', 'bold');
disp('Animación finalizada.');