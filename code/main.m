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

figure("Name", "Sin Refrigeración", "NumberTitle", "off");
plot(t_sr/3600, T_p_sr, "LineWidth", 2, "Color", "b");
title("Temperatura del Panel, sin uso de refrigeración");
xlim([0, 24]); xticks(0:2:24);
grid on;

%% ========== LAZO ABIERTO (ESCALÓN)  ==========

voltaje_bomba_escalon = @(t) 24 * ((9 <= t/3600) & (t/3600 <= 19));
voltaje_ventilador_escalon = @(t) 12 * ((9 <= t/3600) & (t/3600 <= 19));

odefun_la = @(t, y) modelo(t, y, p, ...
    @(t) pert.temperatura_ambiente_despejado(t), ...
    @(t) pert.irradiancia_solar_despejado(t), ...
    @(t) pert.velocidad_viento(t), ...
    0, ...
    voltaje_bomba_escalon, ...
    voltaje_ventilador_escalon, ...
    0);

[t_la, y_la] = ode45(odefun_la, tspan, y0);
T_p_la = y_la(:, 1) - 273.15;

figure("Name", "Lazo Abierto", "NumberTitle", "off");
plot(t_la/3600, T_p_la, "LineWidth", 2, "Color", "b");
title("Temperatura del Panel: Control en Lazo Abierto (Escalón)");
xlim([0, 24]); xticks(0:2:24);
ylim([0, 60]);
grid on;

%% ========== STEP TEST (OBJETIVO 5) - CRITERIO "PRIMER TOQUE" ==========
p_test = p;
Kp_val = p_test.K_p(1);
% Parametros constantes para la prueba
in_amb = @(t) 30; in_irr = @(t) 800; in_wind = @(t) 5.0;

% --- FASE 1: Pre-estabilización (40°C) ---5
p_test.ref = 40 + 273.15;
[t1, y1] = ode45(@(t,y) modelo(t,y,p_test,in_amb,in_irr,in_wind,1,@(t)0,@(t)0,Kp_val), [0 3600], [303.15, 303.15]);

% --- FASE 2: Escalón (55°C) ---
p_test.ref = 55 + 273.15;
[t2, y2] = ode45(@(t,y) modelo(t,y,p_test,in_amb,in_irr,in_wind,1,@(t)0,@(t)0,Kp_val), [3600 10800], y1(end,:));

% Unir datos
t_total = [t1; t2];
T_total = [y1(:,1); y2(:,1)] - 273.15;
T_fase2 = y2(:,1) - 273.15;
t_rel = t2 - t2(1);

% --- CÁLCULO DE MÉTRICAS (MODIFICADO A TU CRITERIO) ---
T_final = mean(T_fase2(end-50:end)); % Valor real donde se estaciona
Delta_V = T_final - T_fase2(1);
banda = 0.01 * abs(Delta_V); % Banda muy estrecha (1%) para detectar el "toque" exacto

% LÓGICA CAMBIADA: Buscamos el PRIMER momento ('first') que entra en la banda
% Esto marca el instante exacto que toca el valor final antes de oscilar
idx_toque = find(abs(T_fase2 - T_final) <= banda, 1, 'first');

if isempty(idx_toque), Ts = 0; else, Ts = t_rel(idx_toque); end

% Otras métricas
ESS = 55 - T_final;
Sobrepaso = max(T_fase2) - T_final;
t_10 = t_rel(find(T_fase2 >= T_fase2(1) + 0.1*Delta_V, 1));
t_90 = t_rel(find(T_fase2 >= T_fase2(1) + 0.9*Delta_V, 1));
Tr = t_90 - t_10;

% --- IMPRIMIR Y GRAFICAR ---
fprintf('\n--- RESULTADOS (Criterio: Primer Toque) ---\n');
fprintf('Estable en: %.2f°C (ESS: %.2f°C)\n', T_final, ESS);
fprintf('Tiempo de Llegada ("Estabilización"): %.1f s\n', Ts);

figure("Name", "Step Test - Primer Toque", "NumberTitle", "off");
% Grafico Temperatura
subplot(1,2,1); plot(t_total/3600, T_total, 'b', 'LineWidth', 2); grid on; hold on;
yline(55, 'r--', 'Ref 55'); yline(T_final, 'k-.', 'Final Real');
xline((3600+Ts)/3600, 'g-', 'Llegada'); % Línea verde marca el momento exacto
title(sprintf('Tiempo hasta tocar valor final: %.1fs', Ts));
xlim([0 2]); ylim([35 60]); xlabel('Tiempo (h)'); ylabel('Temperatura Panel (°C)');

% Grafico Voltaje (Calculado al vuelo para ahorrar código)
u = Kp_val * ((t_total>=3600).*(55+273.15) + (t_total<3600).*(40+273.15) - (T_total+273.15)) + p.offset;
subplot(1,2,2); plot(t_total/3600, min(12, max(0, u)), 'b'); grid on;
title('Voltaje Ventilador');
xlim([0 2]); xlabel('Tiempo (h)');
ylabel('Voltaje (V)'); ylim([-1 13]);



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

for i=1:size(perfil_pert, 1)
    temperatura_ambiente = perfil_pert{i, 1};
    irradiancia_solar = perfil_pert{i, 2};
    velocidad_viento = @(t) pert.velocidad_viento(t);

    nombre_graf = nombre_perfiles(i);

    graficos(nombre_graf, p, temperatura_ambiente, irradiancia_solar, ...
        velocidad_viento, 1, @(t) 0, @(t) 0, Kp_optimo);
end


%% ========== OBJETIVO 6: ESTABILIDAD EN ESCENARIO REAL (DÍA DESPEJADO) ==========
figure("Name", "Comparativa Lazo Cerrado - Día Despejado", "NumberTitle", "off");

% --- Subplot 1: Temperaturas ---
ax1 = subplot(1, 2, 1); hold on; grid on;
title("Respuesta de Temperatura ante diferentes Ganancias (Kp)");
ylabel("Temperatura Panel (°C)");
xlabel("Hora del día");
xlim([0 24]); ylim([0 60]);
yline(55, 'k--', 'Límite 55°C', 'LabelHorizontalAlignment', 'left');

% --- Subplot 2: Esfuerzo de Control ---
ax2 = subplot(1, 2, 2); hold on; grid on;
title("Esfuerzo de Control (Voltaje Ventilador)");
ylabel("Voltaje (V)");
xlabel("Hora del día");
xlim([0 24]); ylim([-1 13]);

% Definición de perturbaciones
f_amb = @(t) pert.temperatura_ambiente_despejado(t);
f_irr = @(t) pert.irradiancia_solar_despejado(t);
f_wind = @(t) pert.velocidad_viento(t);

colores = lines(length(p.K_p));
leyendas = strings(1, length(p.K_p));

% === ARRAY PARA GUARDAR LOS OBJETOS DE LA GRÁFICA ===
graficas = gobjects(1, length(p.K_p));

for i = 1:length(p.K_p)
    Kp_actual = p.K_p(i);

    [t_sim, y_sim] = ode45(@(t,y) modelo(t,y,p, f_amb, f_irr, f_wind, 1, @(t)0, @(t)0, Kp_actual), tspan, y0);

    T_panel = y_sim(:, 1) - 273.15;
    u_calc = Kp_actual * (p.ref - y_sim(:,1)) + p.offset;
    V_vent_plot = min(p.V_vent_MAX, max(u_calc, p.V_MIN));

    % --- Graficamos Temp y GUARDAMOS EL IDENTIFICADOR en 'graficas(i)' ---
    graficas(i) = plot(ax1, t_sim/3600, T_panel, 'LineWidth', 2, 'Color', colores(i,:));

    % --- Graficamos Voltaje ---
    plot(ax2, t_sim/3600, V_vent_plot, 'LineWidth', 1.5, 'Color', colores(i,:));

    leyendas(i) = sprintf('Kp = %.1f', Kp_actual);
end