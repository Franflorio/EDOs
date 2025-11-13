function [x, y, e_local, e_global] = rk4_pvi_auditoria(f, y_exact, x0, y0, h, xfinal)
% rk4_pvi_auditoria  Runge Kutta clasico de orden 4 con auditoria.
%
%   [x, y, e_local, e_global] = rk4_pvi_auditoria(f, y_exact, x0, y0, h, xfinal)

    Nreal = (xfinal - x0) / h;
    N = round(Nreal);

    if abs(Nreal - N) > 1.0e-10
        warning('rk4_pvi_auditoria: Intervalo no multiplo de h. N = %d pasos.', N);
    end

    x = zeros(N + 1, 1);
    y = zeros(N + 1, 1);

    x(1) = x0;
    y(1) = y0;

    e_local = [];
    e_global = [];

    usa_exacta = isa(y_exact, 'function_handle');

    if usa_exacta
        e_local = zeros(N, 1);
        e_global = zeros(N + 1, 1);
        e_global(1) = y_exact(x0) - y0;
    end

    for n = 1:N
        % paso numerico RK4
        k1 = f(x(n), y(n));

        x_med = x(n) + h * 0.5;
        y_med1 = y(n) + h * 0.5 * k1;
        k2 = f(x_med, y_med1);

        y_med2 = y(n) + h * 0.5 * k2;
        k3 = f(x_med, y_med2);

        x_fin = x(n) + h;
        y_fin = y(n) + h * k3;
        k4 = f(x_fin, y_fin);

        x(n + 1) = x_fin;
        y(n + 1) = y(n) + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

        if usa_exacta
            % paso local partiendo de la solucion exacta en x(n)
            y_ex_n = y_exact(x(n));
            y_ex_n1 = y_exact(x(n + 1));

            k1_ex = f(x(n), y_ex_n);

            x_med_ex = x(n) + h * 0.5;
            y_med1_ex = y_ex_n + h * 0.5 * k1_ex;
            k2_ex = f(x_med_ex, y_med1_ex);

            y_med2_ex = y_ex_n + h * 0.5 * k2_ex;
            k3_ex = f(x_med_ex, y_med2_ex);

            x_fin_ex = x(n) + h;
            y_fin_ex = y_ex_n + h * k3_ex;
            k4_ex = f(x_fin_ex, y_fin_ex);

            y_loc_step = y_ex_n + (h / 6) * (k1_ex + 2 * k2_ex + 2 * k3_ex + k4_ex);

            e_local(n) = y_ex_n1 - y_loc_step;
            e_global(n + 1) = y_ex_n1 - y(n + 1);
        end
    end
end

