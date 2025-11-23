%% ========== MODELO (EDO) ==========
function dydt = modelo(t, y, p, temperatura_ambiente, ...
    irradiancia_solar, velocidad_viento, ...
    flag_lazo, v_bomb_fun, v_vent_fun, ...
    K_p)

% Variables de estado en [K]
T_p = y(1);
T_w_in = y(2);

% ========== Entrada ==========

% ---------- LAZO ABIERTO ----------
if flag_lazo == 0
    V_bomb = v_bomb_fun(t);
    V_vent = v_vent_fun(t);


    % ---------- LAZO CERRADO ----------
elseif flag_lazo == 1
    error = p.ref - T_p;
    u_sin_limite = K_p * error + p.offset;

    V_vent = min(p.V_vent_MAX, max(u_sin_limite, p.V_MIN));
    V_bomb = min(p.V_bomb_suficiente, max(u_sin_limite, p.V_MIN));

end

% Variables Externas
v_vient = velocidad_viento(t);
T_amb = temperatura_ambiente(t) + 273.15;
G = irradiancia_solar(t);

% Temperatura del cielo
T_cielo = T_amb - 8;

% Flujo másico
m_dot_w = p.rho_w * p.k_Qf * V_bomb;

% Convección panel solar
if v_vient<5
    h_p = 5.7+3.8*v_vient;
else
    h_p = 6.47+v_vient^(0.78);
end


if V_bomb > 0
    % ---------- BOMBA ENCENDIDA ----------

    T_w_out = (2*m_dot_w * p.C_p_w * T_w_in + p.h_conv_pw * ...
        p.A_c * (2*T_p - T_w_in)) / (2*m_dot_w * p.C_p_w ...
        + p.h_conv_pw * p.A_c);

    % DeltaT_ml sin indeterminaciones
    % Para evitar posibles problemas numéricos
    if abs(T_w_out - T_w_in) > 0 && (T_w_out - T_amb) > 0 && ...
            (T_w_in - T_amb) > 0

        arg_log = (T_w_out - T_amb) / (T_w_in - T_amb);

        if arg_log > 0 && abs(arg_log - 1) > 0
            Delta_T_ml = (T_w_out - T_w_in) / log(arg_log);

        else
            Delta_T_ml = T_w_in - T_amb;
        end
    else
        Delta_T_ml = T_w_in - T_amb;
    end

else
    % ---------- BOMBA APAGADA ----------
    T_w_out = T_w_in;  % Sin flujo, no hay gradiente
    Delta_T_ml = T_w_in - T_amb;  % Para agua estancada, diferencia simple
end


% Flujos de calor
Q_sol = p.A_p * G * (p.alpha_p - p.eta_ref * ...
    (1 + (p.mu_Voc/p.V_mp)*(T_p - p.T_ref)));

Q_conv_p = h_p * p.A_p * (T_p - T_amb);

Q_rad_p = p.epsilon_p * p.sigma * p.A_p * (T_p^4 - T_cielo^4);

Q_cool = m_dot_w * p.C_p_w * (T_w_out - T_w_in);

Q_int = (p.U_base + p.k_vent * V_vent) * p.A_int * Delta_T_ml;

% Ecuaciones diferenciales
dT_p_dt = (Q_sol - Q_conv_p - Q_rad_p - Q_cool) / (p.m_p * p.C_p_p);
dT_w_in_dt = (Q_cool - Q_int) / (p.m_w * p.C_p_w);

% Sistema
dydt = [dT_p_dt;
    dT_w_in_dt];
end