%% ========== PARÁMETROS ==========
function p = parametros()
p = struct();

% Panel Solar
p.m_p = 18.2;           % [kg]
p.A_p = 1.65;           % [m^2]
p.C_p_p = 624;          % [J/kg.K]
p.alpha_p = 0.9;        % Absortividad
p.epsilon_p = 0.88;     % Emisividad

% Parámetros Fotovoltaicos
p.eta_ref = 0.1619;     % Eficiencia referencia
p.mu_Voc = -0.0031;     % [1/K]
p.V_mp = 30.6;          % [V]
p.T_ref = 25 + 273.15;  % [K]

% Sistema de Agua
p.m_w = 5;              % [kg]
p.C_p_w = 4186;         % [J/kg.K]
p.rho_w = 980;          % [kg/m^3]

% Bomba
p.k_Qf = 1.35415e-5;    % [m^3/s.V]

% Interfaz Panel-Agua
p.h_conv_pw = 200.0;    % [W/m^2.K]
p.A_c = 1.3;            % [m^2]

% Intercambiador
p.A_int = 0.78;         % [m^2]
p.U_base = 15.0;        % [W/m^2.K]
p.k_vent = 3.5;         % [W/m^2K.V]

p.sigma = 5.67e-8;


% Parámetros del controlador proporcional
p.K_p = [-0.15, -0.5, -1.0, -2.0];  % Ganancias proporcionales a probar

p.ref = 45 + 273.15;
p.offset = 5;

p.V_MIN = 0;
p.V_vent_MAX = 12;
p.V_bomb_suficiente = 3;

end