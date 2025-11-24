%% ========== GRÁFICOS ==========

function graficos(nombre, p, temperatura_ambiente, irradiancia_solar, ...
    velocidad_viento, flag_lazo, v_vent_fun, v_bomb_fun, ...
    K_p)

% Se soluciona la EDO con las entradas dadas
pert = Entradas();
tspan = 0:60:24*3600;
T0 = pert.temperatura_ambiente_despejado(0) + 273.15;
y0 = [T0, T0];

odefun = @(t,y) modelo(t, y, p, ...
    temperatura_ambiente, ...
    irradiancia_solar, ...
    velocidad_viento,...
    flag_lazo, ...
    v_vent_fun, ...
    v_bomb_fun, ...
    K_p);

[t, y] = ode45(odefun, tspan, y0);
T_p = y(:, 1) - 273.15;

fprintf("Escenario %s: %.2f\n", nombre, max(T_p));

% Para graficar perturbaciones
plot_temp_amb = arrayfun(@(t) temperatura_ambiente(t), t);
plot_irr_sol = arrayfun(@(t) irradiancia_solar(t), t);

% Para graficar entradas
u = K_p * (p.ref - y(:,1)) + p.offset;
salida_vent = min(p.V_vent_MAX, max(u, p.V_MIN));
salida_bomb = min(p.V_bomb_suficiente, max(u, p.V_MIN));

% Máx temperatura alcanzada.
max_temp = max(T_p);

% Para graficar energía consumida
W_vent = salida_vent.^2 / p.R_vent;
W_bomb = salida_bomb.^2 / p.R_bomb;
% Energía consumida en kWh
E_vent = cumtrapz(t, W_vent) / 3.6e6;
E_bomb = cumtrapz(t, W_bomb) / 3.6e6;
E_total = E_vent + E_bomb;

% ---------- Temperatura Panel Solar ----------
figure("Name", compose("Temp/Entradas: %s",nombre), "NumberTitle", "off");
subplot(1, 2, 1);
plot(t/3600, T_p);
xlim([0,24]); xticks(0:2:24);
ylim([0, 60]);
xlabel("Hora del día (hrs)"); ylabel("Temperatura del Panel (°C)");
yline(max_temp, 'r--', compose("Máx: %.2f°C", max_temp),'Color', 'r', 'LineWidth', 2);
grid on;

% ---------- Voltajes Ventilador/Bomba ----------
subplot(1, 2, 2);
plot(t/3600, salida_vent, "LineStyle", "--", "LineWidth", 2, "Color", "b");
hold on; grid on;
plot(t/3600, salida_bomb, "LineStyle", "-.", "LineWidth", 2, "Color", "r");
legend("V. Ventilador", "V. Bomba", "Location", "northwest");
xlim([0,24]); xticks(0:2:24);
ylim([-1, 13]);
xlabel("Hora del día (hrs)"); ylabel("Voltaje (V)");
hold off;

% ---------- Energía Consumida ----------
figure("Name", compose("Energía Consumida: %s", nombre), "NumberTitle", "off");
plot(t/3600, E_vent, "LineWidth", 2, "Color", "b");
hold on; grid on;
plot(t/3600, E_bomb, "LineWidth", 2, "Color", "r");
plot(t/3600, E_total, "LineWidth", 2, "Color", "g");
legend("E. consum. Ventilador", "E. consum. Bomba", "E. Total consum.", "Location", "northwest");
xlim([0,24]); xticks(0:2:24);
ylim([-0.001, 0.04]); yticks(0:0.01:0.04);
xlabel("Hora del día (hrs)"); ylabel("Energía Consumida (kWh)");

% ---------- Temperatura Ambiente ----------
figure("Name", compose("Perturbaciones: %s", nombre), "NumberTitle", "off");
subplot(1, 2, 1);
plot(t/3600, plot_temp_amb, "LineWidth", 1, "Color", "g");
xlim([0,24]); xticks(0:2:24);
ylim([0, 40]);
xlabel("Hora del día (hrs)"); ylabel("Temperatura del Panel (°C)");
grid on;

% ---------- Irradiancia Solar ----------
subplot(1, 2, 2);
plot(t/3600, plot_irr_sol, "LineWidth", 1, "Color", "y");
xlim([0,24]); xticks(0:2:24);
ylim([-100, 1100]);
xlabel("Hora del día (hrs)"); ylabel("Irradiancia (W/m^2)");
grid on;

end