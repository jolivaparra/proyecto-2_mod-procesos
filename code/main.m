%% ========== CÓDIGO PRINCIPAL ==========

clear; close all; clc;

% =========== Carga de Parámetros y Perturbaciones ==========
p = parametros();
pert = Entradas();

tspan = 0:69:24*3600;
T0 = pert.temperatura_ambiente_despejado(0) + 273.15;
y0 = [T0, T0];

%% ========== SIN REFRIGERACION ==========

odefun_snrefr = @(t, y) modelo(t, y, p, ...
    @(t) pert.temperatura_ambiente_despejado(t), ...
    @(t) pert.irradiancia_solar_despejado(t), ...
    @(t) pert.velocidad_viento(t), ...
    0, ...
    @(t) 0, ...
    @(t) 0, ...
    0);

[t_sr, y_sr] = ode45(odefun_snrefr, tspan, y0);
T_p_sr = y_sr(:, 1) - 273.15;

fprintf("Máx día desjepado sin refrigeración: %.2f\n", max(T_p_sr));

figure("Name", "Sin Refrigeración", "NumberTitle", "off");
plot(t_sr/3600, T_p_sr, "LineWidth", 2, "Color", "b");
xlabel("Hora del día (hrs)"); ylabel("Temperatura del Panel (°C)");
xlim([0, 24]); xticks(0:2:24);
grid on;


%% ========== LAZO ABIERTO (ESCALÓN)  ==========

voltaje_bomba_escalon = @(t) 24 * ((9 <= t/3600) & (t/3600 <= 19));
voltaje_ventilador_escalon = @(t) 12 * ((9 <= t/3600) & (t/3600 <= 19));

odefun_la = @(t, y) modelo(t, y, p, ...
    @(t) pert.temperatura_ambiente_despejado(t), ...
    @(t) pert.irradiancia_solar_despejado(t), ...
    @(t) pert.velocidad_viento(t), 0, voltaje_bomba_escalon, ...
    voltaje_ventilador_escalon, 0);

[t_la, y_la] = ode45(odefun_la, tspan, y0);
T_p_la = y_la(:, 1) - 273.15;

fprintf("Máx día desjepado con entrada escalón (Lazo Abierto): %.2f\n", max(T_p_la));

figure("Name", "Lazo Abierto", "NumberTitle", "off");
plot(t_la/3600, T_p_la, "LineWidth", 2, "Color", "b");
xlabel("Hora del día (hrs)"); ylabel("Temperatura del Panel (°C)");
xlim([0, 24]); xticks(0:2:24);
ylim([0, 60]);
grid on;

%% ========== OBJETIVO 6: ESTABILIDAD EN ESCENARIO REAL (DÍA DESPEJADO) ==========
figure("Name", "Comparativa Kp - Día Despejado", "NumberTitle", "off");

% --- Subplot 1: Temperaturas ---
ax1 = subplot(1, 2, 1); hold on; grid on;
ylabel("Temperatura Panel (°C)");
xlabel("Hora del día (hrs)");
xlim([0 24]); ylim([0 60]);
yline(55, 'k--', 'Límite 55°C', 'LabelHorizontalAlignment', 'right');

% --- Subplot 2: Esfuerzo de Control ---
ax2 = subplot(1, 2, 2); hold on; grid on;
ylabel("Voltaje (V)");
xlabel("Hora del día (hrs)");
xlim([0 24]); ylim([-1 13]);

% Definición de perturbaciones
f_amb = @(t) pert.temperatura_ambiente_despejado(t);
f_irr = @(t) pert.irradiancia_solar_despejado(t);
f_wind = @(t) pert.velocidad_viento(t);

colores = lines(length(p.K_p));
leyendas = strings(1, length(p.K_p));

% === ARRAYS PARA GUARDAR LOS OBJETOS DE LA GRÁFICA ===
graficas_temp = gobjects(1, length(p.K_p));
graficas_volt = gobjects(1, length(p.K_p));

for i = 1:length(p.K_p)
    Kp_actual = p.K_p(i);

    [t_sim, y_sim] = ode45(@(t,y) modelo(t,y,p, f_amb, f_irr, f_wind, 1, @(t)0, @(t)0, Kp_actual), tspan, y0);

    T_panel = y_sim(:, 1) - 273.15;
    u_calc = Kp_actual * (p.ref - y_sim(:,1)) + p.offset;
    V_vent_plot = min(p.V_vent_MAX, max(u_calc, p.V_MIN));

    % --- Graficamos Temp y GUARDAMOS EL OBJETO ---
    graficas_temp(i) = plot(ax1, t_sim/3600, T_panel, 'LineWidth', 2, 'Color', colores(i,:));

    % --- Graficamos Voltaje y GUARDAMOS EL OBJETO ---
    % <--- MODIFICADO: Ahora guardamos esto en 'graficas_volt(i)'
    graficas_volt(i) = plot(ax2, t_sim/3600, V_vent_plot, 'LineWidth', 1.5, 'Color', colores(i,:));

    % Creamos el texto de la leyenda para este Kp
    leyendas(i) = compose('K_p = %.2f', Kp_actual);
end

% === AÑADIR LAS LEGENDS AL FINAL ===
% <--- AGREGADO: Aquí aplicamos las leyendas a los objetos guardados
legend(ax1, graficas_temp, leyendas, 'Location', 'best', 'Interpreter', 'tex');
legend(ax2, graficas_volt, leyendas, 'Location', 'northeast', 'Interpreter', 'tex');

%% ========== STEP TEST (OBJETIVO 5) - FINAL ==========
p_test = p;
Kp_val = p_test.K_p(1);
% Parametros constantes
in_amb = @(t) 30; in_irr = @(t) 800; in_wind = @(t) 5.0;

% --- FASE 1: Pre-estabilización (1 hora a 40°C) ---
p_test.ref = 40 + 273.15;
[t1, y1] = ode45(@(t,y) modelo(t,y,p_test,in_amb,in_irr,in_wind,1,@(t)0,@(t)0,Kp_val), ...
    [0 3600], [303.15, 303.15]);

% --- FASE 2: Escalón (1 hora a 55°C) ---
p_test.ref = 55 + 273.15;
[t2, y2] = ode45(@(t,y) modelo(t,y,p_test,in_amb,in_irr,in_wind,1,@(t)0,@(t)0,Kp_val), ...
    [3600 7200], y1(end,:));

% Unir datos
t_total = [t1; t2];
T_total = [y1(:,1); y2(:,1)] - 273.15;
T_fase2 = y2(:,1) - 273.15;
t_rel = t2 - t2(1);

% --- CÁLCULO DE MÉTRICAS ---

% 1. Centro de Oscilación (Último 30% de datos)
idx_cola = round(length(T_fase2) * 0.7) : length(T_fase2);
datos_cola = T_fase2(idx_cola);
T_final_real = (max(datos_cola) + min(datos_cola)) / 2;

% 2. Delta del viaje
Delta_V = T_final_real - T_fase2(1);

% 3. Banda de Tolerancia (5%)
banda = 0.05 * abs(Delta_V);

% 4. Tiempo de Estabilización (Ts) - Primer Toque
idx_toque = find(abs(T_fase2 - T_final_real) <= banda, 1, 'first');
if isempty(idx_toque), Ts = 0; else, Ts = t_rel(idx_toque); end

% 5. Tiempo de Subida (Tr) - 10% a 90%
val_10 = T_fase2(1) + 0.1 * Delta_V;
val_90 = T_fase2(1) + 0.9 * Delta_V;
idx_10 = find(T_fase2 >= val_10, 1, 'first');
idx_90 = find(T_fase2 >= val_90, 1, 'first');
if ~isempty(idx_10) && ~isempty(idx_90)
    Tr = t_rel(idx_90) - t_rel(idx_10);
else
    Tr = 0;
end

% 6. Otras
ESS = 55 - T_final_real;
Sobrepaso = max(T_fase2) - T_final_real;

% --- IMPRIMIR RESULTADOS ---
fprintf('\n--- RESULTADOS STEP TEST ---\n');
fprintf('Centro Oscilación: %.4f°C\n', T_final_real);
fprintf('Banda (5%%): +/- %.4f°C\n', banda);
fprintf('Ts: %.1f s | Tr: %.1f s | ESS: %.2f °C | Sobrepaso: %.2f\n\n', Ts, Tr, ESS, Sobrepaso);

% --- GRÁFICOS FINAL ---
figure("Name", "Step Test - Métricas", "NumberTitle", "off");

% Subplot 1: Temperatura
subplot(1,2,1);
plot(t_total/3600, T_total, 'b', 'LineWidth', 2); grid on; hold on;
yline(55, 'r--', 'Ref 55'); yline(45, 'r--', 'Ref 45');
yline(T_final_real, 'k-.', 'Centro Osc.');

% Bandas oscuras y visibles
yline(T_final_real + banda, 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);
yline(T_final_real - banda, 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);

% Línea de Tiempo de Estabilización
if Ts > 0
    % (3600 + Ts) porque el escalón empieza en la hora 1
    xline((3600+Tr)/3600, 'm-', compose('Tr=%.0fs', Ts), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', 'b');
    xline((3600+Ts)/3600, 'm-', compose('Ts=%.0fs', Ts), 'LineWidth', 2);
end

title("Respuesta Transitoria");
xlim([0 2]);
ylim([40 57.5]);
xlabel('Tiempo (h)'); ylabel('Temperatura del Panel (°C)');

% Subplot 2: Esfuerzo de Control
subplot(1,2,2);
u_calc = Kp_val * ((t_total>=3600).*(55+273.15) + (t_total<3600).*(40+273.15) - (T_total+273.15)) + p.offset;
V_vent_plot = min(12, max(0, u_calc));
plot(t_total/3600, V_vent_plot, 'b', 'LineWidth', 1.5);
grid on;
title('Esfuerzo de Control');
xlim([0 2]); % Límite Fijo solicitado
ylim([-1 13]);
xlabel('Tiempo (h)'); ylabel('Voltaje (V)');

%% ========== LAZO CERRADO (VALIDACIÓN DEL DISEÑO FINAL) ==========
% Objetivo: Demostrar que el Kp elegido (-0.1) funciona en los 3 climas.
% NO comparamos Kp aquí, solo validamos el comportamiento.

nombre_perfiles = ["Día Despejado", "Día Nublado", "Día Intermitente"];

% Array de funciones de pertubación
perfil_pert = {
    @(t) pert.temperatura_ambiente_despejado(t), @(t) pert.irradiancia_solar_despejado(t);
    @(t) pert.temperatura_ambiente_nublado(t),   @(t) pert.irradiancia_solar_nublado(t);
    @(t) pert.temperatura_ambiente_intermitente(t), @(t) pert.irradiancia_solar_intermitente(t)
    };

Kp_optimo = p.K_p(1);

fprintf("--- Max Temp con K_p:%.2f ---\n", Kp_optimo);

for i=1:size(perfil_pert, 1)
    temperatura_ambiente = perfil_pert{i, 1};
    irradiancia_solar = perfil_pert{i, 2};
    velocidad_viento = @(t) pert.velocidad_viento(t);

    nombre_graf = nombre_perfiles(i);

    graficos(nombre_graf, p, temperatura_ambiente, irradiancia_solar, ...
        velocidad_viento, 1, @(t) 0, @(t) 0, Kp_optimo);
end


