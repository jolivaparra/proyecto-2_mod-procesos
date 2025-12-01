%% ========== GRÁFICOS ==========
function graficos(nombre, p, temperatura_ambiente, ...
    irradiancia_solar, velocidad_viento, flag_lazo, ...
    v_vent_fun, v_bomb_fun, K_c)

% Soluciona la EDO con las entradas dadas
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
    K_c);

[t, y] = ode45(odefun, tspan, y0);
T_p = y(:, 1) - 273.15;
t_h = t/3600;

% Grafica entradas
u = K_c * (p.ref - y(:,1)) + p.offset;
salida_vent = min(p.V_vent_MAX, max(u, p.V_MIN));
salida_bomb = min(p.V_bomb_suficiente, max(u, p.V_MIN));

% Máx temperatura
max_temp = max(T_p);

% ---------- Energía Consumida ----------
% Calcula energía consumida acomulada para graficar
W_vent = salida_vent.^2 / p.R_vent;
W_bomb = salida_bomb.^2 / p.R_bomb;
% Calcula energía consumida en kWh
E_vent = cumtrapz(t, W_vent) / 3.6e+6;
E_bomb = cumtrapz(t, W_bomb) / 3.6e+6;
E_total = E_vent + E_bomb;

% Imprime en la consola datos térmicos y energéticos relevantes
fprintf("%s: \t %.2f °C | %.3f kWh | %.3f kWh | %.3f kWh \n", ...
    nombre, max_temp, E_vent(end), E_bomb(end), E_total(end));

% Crea vector para graficar perturbaciones
temp_amb = arrayfun(@(t) temperatura_ambiente(t), t);
irr_sol = arrayfun(@(t) irradiancia_solar(t), t);

% ---------- Temperatura Panel Solar ----------
figure("Name", compose("Tp/Entradas: %s",nombre), "NumberTitle", "off");
eje1 = subplot(1, 2, 1);
plot(t_h, T_p);
yline(max_temp, 'r--', compose("Máx: %.2f°C", max_temp),'Color', 'r', 'LineWidth', 2);
ylim([0, 60]); ylabel("Temperatura del Panel (°C)");
grid on;

% ---------- Voltajes Ventilador/Bomba ----------
eje2 = subplot(1, 2, 2);
plot(t_h, salida_vent, "LineStyle", "--", "LineWidth", 2, "Color", "b");
hold on; grid on;
plot(t_h, salida_bomb, "LineStyle", "-.", "LineWidth", 2, "Color", "r");
legend("V. Ventilador", "V. Bomba", "Location", "northwest");
ylim([-1, 13]); ylabel("Voltaje (V)");
hold off;

if nombre=="Despejado"
    % ---------- Energía Consumida ----------
    figure("Name", compose("E. Consum.: %s", nombre), "NumberTitle", "off");
    eje3 = gca;
    plot(t_h, E_vent, "LineWidth", 2, "Color", "b");
    hold on; grid on;
    plot(t_h, E_bomb, "LineWidth", 2, "Color", "r");
    plot(t_h, E_total, "LineWidth", 2, "Color", "g");
    legend("E. consum. Ventilador", ...
        "E. consum. Bomba", ...
        "E. Total consum.", ...
        "Location", "northwest");
    ylim([-0.001, 0.04]); yticks(0:0.01:0.04);
    ylabel("Energía Consumida (kWh)");
end

% ---------- Temperatura Ambiente ----------
figure("Name", compose("Pert: %s", nombre), "NumberTitle", "off");
eje4 = subplot(1, 2, 1);
plot(t_h, temp_amb, "LineWidth", 2, "Color", [0,0.5,0]);
ylim([0, 40]); ylabel("Temperatura del Panel (°C)");
grid on;

% ---------- Irradiancia Solar ----------
eje5 = subplot(1, 2, 2);
plot(t_h, irr_sol, "LineWidth", 1.75, "Color", [0.85,0.85,0]);
ylim([-100, 1100]); ylabel("Irradiancia (W/m^2)");
grid on;

% === CONFIGURACIONES GENERALES DE LOS GRÁFICOS ===
% Agrupa todos los gráficos
ejes = [eje1, eje2, eje4, eje5];
if nombre=="Despejado"
    ejes = [ejes, eje3];
end

% Añade datos comunes al eje X
set(ejes, ...
    'xlim', [0, 24], ...
    'xtick', 0:2:24);
xlabel(ejes, "Hora del día (hrs)")
end