%% ========== ENTRADAS DE PERTURBACIONES ==========z
classdef Entradas
    methods

        % Velocidad del viento (Misma para todos los perfiles)
        function v_vient = velocidad_viento(~, t)
            v_vient = 5+3*sin(2*pi*t/3600*24);
        end

        % ========== DÍA DESPEJADO ==========

        % Irradiancia solar
        function G = irradiancia_solar_despejado(~, t)

            G_max = 1000;      % Irradiancia pico (W/m^2)
            t_amanecer = 6.5;    % Hora del amanecer (horas)
            t_atardecer = 20.5;  % Hora del ocaso (horas)
            t_h = t / 3600;

            arg = pi * (t_h-t_amanecer)/(t_atardecer-t_amanecer);

            G = (t_amanecer<=t_h & t_h<=t_atardecer).*G_max.*sin(arg);

        end

        % Temperatura ambiente
        function T_amb = temperatura_ambiente_despejado(~, t)

            t_h = t/3600;

            t_min = 6.5;
            t_sunset = 20.5;

            Y = t_sunset-t_min;
            Z = 24-Y;

            T_x = 35.2;
            T_N = 10;
            T_s = 25.1018;

            a = 1.8;
            b = 2.5;

            if t_h>=t_min && t_h<=t_sunset

                T_amb = (T_x-T_N)*sin(pi*(t_h-t_min)/(Y+2*a))+T_N;

            else

                % Ajuste de continuidad y periodicidad
                % para comenzar simulaciones en tiempo 00:00
                if t_h<t_min
                    T_amb = T_N+(T_s-T_N)*exp(-(b*((t_h+24)-t_sunset))/(Z));

                else
                    T_amb = T_N+(T_s-T_N)*exp(-(b*(t_h-t_sunset))/(Z));
                end
            end
        end

        % ========== DÍA NUBLADO ==========

        % Irradiancia solar
        function G = irradiancia_solar_nublado(~, t)

            G_max = 1000*0.5;      % Irradiancia pico (W/m^2)
            t_amanecer = 6.5;    % Hora del amanecer (horas)
            t_atardecer = 20.5;  % Hora del ocaso (horas)
            t_h = t / 3600;

            arg = pi * (t_h-t_amanecer)/(t_atardecer-t_amanecer);

            G = (t_amanecer<=t_h & t_h<=t_atardecer).*G_max.*sin(arg);

        end

        % Temperatura ambiente
        function T_amb = temperatura_ambiente_nublado(~, t)

            t_h = t/3600;

            t_min = 6.5;
            t_sunset = 20.5;

            Y = t_sunset-t_min;
            Z = 24-Y;

            T_x = 25;
            T_N = 12;      % <-- Mínima similar o más alta
            T_s = 20;      % <-- Temp. atardecer

            a = 1.8;
            b = 2.5;

            if t_h>=t_min && t_h<=t_sunset

                T_amb = (T_x-T_N)*sin(pi*(t_h-t_min)/(Y+2*a))+T_N;

            else

                % Ajuste de continuidad y periodicidad
                % para comenzar simulaciones en tiempo 00:00
                if t_h<t_min
                    T_amb = T_N+(T_s-T_N)*exp(-(b*((t_h+24)-t_sunset))/(Z));

                else
                    T_amb = T_N+(T_s-T_N)*exp(-(b*(t_h-t_sunset))/(Z));
                end
            end
        end

        % ========== DÍA PARCIALMENTE NUBLADO ==========
        function G = irradiancia_solar_intermitente(obj, t)

            % Base
            G_despejado = obj.irradiancia_solar_despejado(t);

            % Parámetros de Nubes (Locales)
            periodo_1 = 2.5 * 3600;  % Nube grande cada 2.5 horas
            duracion_1 = 25 * 60;    % Dura 25 minutos
            atenuacion_1 = 0.3;      % Pasa el 30%

            periodo_2 = 1.1 * 3600;  % Nube mediana cada 1.1 horas
            duracion_2 = 8 * 60;     % Dura 8 minutos
            atenuacion_2 = 0.6;      % Pasa el 60%

            % Lógica de Pulso Abrupta (rem)
            if rem(t, periodo_1) < duracion_1
                factor_1 = atenuacion_1;
            else
                factor_1 = 1.0;
            end

            desfase_2 = 15 * 60;
            if rem(t + desfase_2, periodo_2) < duracion_2
                factor_2 = atenuacion_2;
            else
                factor_2 = 1.0;
            end

            factor_nube_total = factor_1 * factor_2;
            G = G_despejado * factor_nube_total;
        end

        % 2. TEMPERATURA AMBIENTE (CON INERCIA Y CORRELACIÓN SOLAR)
        %    Reacciona LENTAMENTE solo al clima general (nube grande)
        %    y SOLO si el sol está fuera.

        function T_amb = temperatura_ambiente_intermitente(obj, t)

            % 1. Obtenemos la temperatura base (día despejado)
            T_despejado = obj.temperatura_ambiente_despejado(t);

            % 2. Verificamos si es de día
            %    (Usamos la irradiancia base como interruptor)
            G_despejado = obj.irradiancia_solar_despejado(t);

            if G_despejado <= 0
                % Si es de noche, NO hay efecto de nube.
                reduccion = 0;

            else
                % Si es de día, aplicamos la inercia térmica.

                % --- Parámetros de Inercia (Locales) ---

                % El "clima" general sigue el ciclo de nubes más largo (2.5h)
                periodo_clima = 2.5 * 3600;

                % La T_amb bajará suavemente un máximo de 4 grados
                max_reduccion_C = 4.0;

                % Desfase por Inercia (45 min de retraso)
                retraso = 0.75 * 3600;

                % --- Cálculo de la Onda Suave ---

                % Onda cos() en la misma frecuencia que la nube grande
                onda_clima = cos( (2 * pi / periodo_clima) * (t - retraso) );

                % Mapeamos [-1, 1] -> [0, 1]
                factor_inercial_suave = (onda_clima + 1) / 2;

                % Invertimos: [1 (frío), 0 (cálido)]
                reduccion_norm = 1 - factor_inercial_suave;

                reduccion = reduccion_norm * max_reduccion_C;
            end

            % 3. Aplicamos la reducción
            T_amb = T_despejado - reduccion;
        end

    end
end