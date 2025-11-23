%% ========== GRÁFICOS ==========

function graficos(nombre, p, temperatura_ambiente, irradiancia_solar, ...
                    velocidad_viento, flag_lazo, v_vent_fun, v_bomb_fun, ...
                    K_p)


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

    plot_temp_amb = arrayfun(@(t) temperatura_ambiente(t), t);
    plot_irr_sol = arrayfun(@(t) irradiancia_solar(t), t);

    salida_sin_limite = K_p * (p.ref - y(:,1)) + p.offset;
    salida_vent = min(p.V_vent_MAX, max(salida_sin_limite, p.V_MIN));
    salida_bomb = min(p.V_bomb_suficiente, max(salida_sin_limite, p.V_MIN));


    % ---------- Temperatura Panel Solar ----------
    figure("Name", nombre, "NumberTitle", "off");
    subplot(2, 2, 1);
    plot(t/3600, T_p);
    title("Temperatura del Panel Solar");
    xlim([0,24]); xticks(0:2:24);
    ylim([0, 60]);
    xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");
    grid on;
    
    % ---------- Voltajes Ventilador/Bomba ----------
    subplot(2, 2, 2);
    plot(t/3600, salida_vent, "LineStyle", "--", "LineWidth", 2, "Color", "b");
    hold on; grid on;
    plot(t/3600, salida_bomb, "LineStyle", "-.", "LineWidth", 2, "Color", "r");
    legend("V. Ventilador", "V. Bomba", "Location", "northwest");
    title("Voltaje aplicado a la bomba y ventilador");
    xlim([0,24]); xticks(0:2:24);
    ylim([-1, 13]);
    xlabel("Tiempo (h)"); ylabel("Voltaje (V)");

    % ---------- Temperatura Ambiente ----------

    subplot(2, 2, 3);
    plot(t/3600, plot_temp_amb, "LineWidth", 1, "Color", "g");
    title("Temperatura Ambiente");
    xlim([0,24]); xticks(0:2:24);
    ylim([0, 40]);
    xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");
    grid on;

    % ---------- Irradiancia Solar ----------
    subplot(2, 2, 4);
    plot(t/3600, plot_irr_sol, "LineWidth", 1, "Color", "y");
    title("Irradiancia Solar");
    xlim([0,24]); xticks(0:2:24);
    ylim([-100, 1100]);
    xlabel("Tiempo (h)"); ylabel("Irradiancia (W/m^2)");
    grid on;

end