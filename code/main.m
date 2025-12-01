%% ========== CÓDIGO PRINCIPAL ==========
clear; close all; clc;

% =========== Carga de Parámetros y Perturbaciones ==========
p = parametros();
pert = Entradas();

tspan = 0:69:24*3600;
T0 = pert.temperatura_ambiente_despejado(0) + 273.15;
y0 = [T0, T0];

%% ========== LAZO ABIERTO (ESCALÓN)  ==========
% Define entradas y resuelve EDO
voltaje_bomba_escalon = @(t) 24 * ((9 <= t/3600) & (t/3600 <= 19));
voltaje_ventilador_escalon = @(t) 12 * ((9 <= t/3600) & (t/3600 <= 19));

odefun_la = @(t, y) modelo(t, y, p, ...
    @(t) pert.temperatura_ambiente_despejado(t), ...
    @(t) pert.irradiancia_solar_despejado(t), ...
    @(t) pert.velocidad_viento(t), 0, voltaje_bomba_escalon, ...
    voltaje_ventilador_escalon, 0);

[t_la, y_la] = ode45(odefun_la, tspan, y0);
T_p_la = y_la(:, 1) - 273.15;
t_la_h = t_la/3600;

% Máx temperatura alcanzada.
max_temp = max(T_p_la);

% ---------- Energía Consumida ----------
% Calcula energía consumida acomulada para graficar
W_vent = arrayfun(voltaje_ventilador_escalon, t_la).^2 / p.R_vent;
W_bomb = arrayfun(voltaje_bomba_escalon, t_la).^2 / p.R_bomb;
% Calcula energía consumida en kWh
E_vent = cumtrapz(t_la, W_vent) / 3.6e+6;
E_bomb = cumtrapz(t_la, W_bomb) / 3.6e+6;
E_total = E_vent + E_bomb;

% Imprime en la consola datos térmicos y energéticos relevantes
fprintf("Lazo Abierto:\t Temp Máx |  E_vent   |   E_bomb  | E_total \n")
fprintf("\t\t %.2f °C | %.3f kWh | %.3f kWh | %.3f kWh \n", max_temp, E_vent(end), E_bomb(end), E_total(end));

% Grafica Temperatura del Panel en Lazo Abierto
figure("Name", "Lazo Abierto", "NumberTitle", "off");
plot(t_la_h, T_p_la, "LineWidth", 2, "Color", "b");
xlabel("Hora del día (hrs)"); ylabel("Temperatura del Panel (°C)");
xlim([0, 24]); xticks(0:2:24);
ylim([0, 60]);
grid on;

% Grafica Energía consumida por actuadores
figure("Name", "E Consum. Lazo Abierto", "NumberTitle", "off");
hold on; grid on;
plot(t_la_h, E_total, "LineWidth", 2, "Color", "g");
plot(t_la_h, E_vent, "LineWidth", 2, "Color", "b");
plot(t_la_h, E_bomb, "LineWidth", 2, "Color", "r");
legend("E. consum. Ventilador", "E. consum. Bomba", "E. Total consum.", "Location", "northwest");
xlim([0,24]); xticks(0:2:24);
ylim([-0.1, 0.8]); yticks(0:0.2:0.8);
xlabel("Hora del día (hrs)"); ylabel("Energía Consumida (kWh)");

%% ========== VARIACIONES DE Kc ==========
% Prepara gráficos
figure("Name", "Comparativa Kc - Día Despejado", "NumberTitle", "off");

% --- Subplot 1: Temperaturas ---
ax1 = subplot(1, 2, 1); hold on; grid on;
ylabel("Temperatura Panel (°C)");
xlabel("Hora del día (hrs)");
xlim([0, 24]); xticks(0:2:24);
ylim([0 60]);

yline(55, 'k--', 'Límite 55°C', 'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

% --- Subplot 2: Esfuerzo de Control ---
ax2 = subplot(1, 2, 2); hold on; grid on;
ylabel("Voltaje (V)");
xlabel("Hora del día (hrs)");
xlim([0, 24]); xticks(0:2:24);
ylim([-1 13]);

% Definición de perturbaciones
f_amb = @(t) pert.temperatura_ambiente_despejado(t);
f_irr = @(t) pert.irradiancia_solar_despejado(t);
f_wind = @(t) pert.velocidad_viento(t);

colores = lines(length(p.K_c));

% --- BUCLE DE SIMULACIÓN ---
for i = 1:length(p.K_c)
    Kc_actual = p.K_c(i);

    % Resuelve EDO's
    [t_sim, y_sim] = ode45(@(t,y) modelo(t,y,p, f_amb, f_irr, f_wind, 1, @(t)0, @(t)0, Kc_actual), tspan, y0);

    % Calcula variables a graficar
    T_panel = y_sim(:, 1) - 273.15;
    u_calc = Kc_actual * (p.ref - y_sim(:,1)) + p.offset;
    V_vent_plot = min(p.V_vent_MAX, max(u_calc, p.V_MIN));

    % Crea leyenda
    etiqueta_leyenda = sprintf('K_c = %.2f', Kc_actual);

    % --- Graficar Temperatura ---
    plot(ax1, t_sim/3600, T_panel, ...
        'LineWidth', 2, ...
        'Color', colores(i,:), ...
        'DisplayName', etiqueta_leyenda);

    % --- Graficar Voltaje ---
    plot(ax2, t_sim/3600, V_vent_plot, ...
        'LineWidth', 1.5, ...
        'Color', colores(i,:), ...
        'DisplayName', etiqueta_leyenda);
end

% === ACTIVA LEYENDAS ===
legend(ax1, 'Location', 'best', 'Interpreter', 'tex');
legend(ax2, 'Location', 'northeast', 'Interpreter', 'tex');

%% ========== USO DE MÉTRICAS (STEP TEST) ==========
p_test = p;
Kc_val = p_test.K_c(1);

% Parametros constantes
in_amb = @(t) 30; in_irr = @(t) 800; in_wind = @(t) 5.0;

% --- FASE 1: Pre-estabilización (1 hora a 40°C) ---
p_test.ref = 40 + 273.15;
[t1, y1] = ode45(@(t,y) modelo(t,y,p_test,in_amb,in_irr,in_wind,1,@(t)0,@(t)0,Kc_val), ...
    [0 3600], [273.15 + 30, 273.15 + 30]);

% --- FASE 2: Escalón (1 hora a 55°C) ---
p_test.ref = 55 + 273.15;
[t2, y2] = ode45(@(t,y) modelo(t,y,p_test,in_amb,in_irr,in_wind,1,@(t)0,@(t)0,Kc_val), ...
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
% Promedio de Osc_max y Osc_min
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

% Error de estado estacionario
ESS = 55 - T_final_real;
% Sobrepaso
Sobrepaso = max(T_fase2) - T_final_real;

% --- IMPRIME RESULTADOS EN CONSOLA ---
fprintf('\n--- RESULTADOS ENTORNO CONTROLADO ---\n');
fprintf('Tr: %.1f s | Ts: %.1f s | ESS: %.2f °C | Sobrepaso: %.2f\n\n', Tr, Ts, ESS, Sobrepaso);
fprintf('Centro Oscilación: %.4f°C\n', T_final_real);
fprintf('Banda (5%%): +/- %.4f°C\n', banda);

% --- GRÁFICOS FINALES ---
figure("Name", "Ent. Controlado: Métricas", "NumberTitle", "off");

% Subplot 1: Temperatura
subplot(1,2,1);
plot(t_total/3600, T_total, 'b', 'LineWidth', 2); grid on; hold on;
xlim([0 2]);
ylim([40 57.5]);
xlabel('Tiempo (h)'); ylabel('Temperatura del Panel (°C)');

% Referencia y valor final
yline(55, 'r--', 'Ref 55'); yline(45, 'r--', 'Ref 45');
yline(T_final_real, 'k-.', 'Centro Osc.');

% Bandas oscuras y visibles
yline(T_final_real + banda, 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);
yline(T_final_real - banda, 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);

% Línea de Tiempo de Estabilización
xline((3600+Tr)/3600, 'm-', compose('Tr=%.0fs', Ts), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', 'b');
xline((3600+Ts)/3600, 'm-', compose('Ts=%.0fs', Ts), 'LineWidth', 2);

% Subplot 2: Esfuerzo de Control
subplot(1,2,2);
u_calc = Kc_val * ((t_total>=3600).*(55+273.15) + (t_total<3600).*(40+273.15) - (T_total+273.15)) + p.offset;
V_vent_plot = min(12, max(0, u_calc));
plot(t_total/3600, V_vent_plot, 'b', 'LineWidth', 1.5);
grid on;
xlim([0 2]); % Límite Fijo solicitado
ylim([-1 13]);
xlabel('Tiempo (h)'); ylabel('Voltaje (V)');

%% ========== LAZO CERRADO: VALIDACIÓN DEL DISEÑO FINAL (Kc = -0.25) ==========

nombre_perfiles = ["Despejado", "Nublado", "Intermitente"];

% Array de funciones de pertubación
perfil_pert = {
    @(t) pert.temperatura_ambiente_despejado(t), @(t) pert.irradiancia_solar_despejado(t);
    @(t) pert.temperatura_ambiente_nublado(t), @(t) pert.irradiancia_solar_nublado(t);
    @(t) pert.temperatura_ambiente_intermitente(t), @(t) pert.irradiancia_solar_intermitente(t)
    };

Kp_optimo = p.K_c(1);

fprintf("--- (K_c:%.3f) Max Temp |  E_vent   |  E_bomb   |  E_total \n", Kp_optimo);

for i=1:size(perfil_pert, 1)
    temperatura_ambiente = perfil_pert{i, 1};
    irradiancia_solar = perfil_pert{i, 2};
    velocidad_viento = @(t) pert.velocidad_viento(t);

    nombre_graf = nombre_perfiles(i);

    % Grafica Temperatura del Panel, Entradas y Perturbaciones Aplicadas
    graficos(nombre_graf, p, temperatura_ambiente, irradiancia_solar, ...
        velocidad_viento, 1, @(t) 0, @(t) 0, Kp_optimo);
end