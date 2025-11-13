function [x, y, e_local, e_global] = euler_mejorado_pvi_auditoria(f, y_exact, x0, y0, h, xfinal)
% euler_mejorado_pvi_auditoria  Euler mejorado (Heun) con auditoria de errores.
%
%   [x, y, e_local, e_global] = euler_mejorado_pvi_auditoria(f, y_exact, x0, y0, h, xfinal)

    Nreal = (xfinal - x0) / h;
    N = round(Nreal);

    if abs(Nreal - N) > 1.0e-10
        warning('euler_mejorado_pvi_auditoria: Intervalo no multiplo de h. N = %d pasos.', N);
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
        % paso numerico
        k1 = f(x(n), y(n));
        x_next = x(n) + h;
        y_pred = y(n) + h * k1;
        k2 = f(x_next, y_pred);

        y(n + 1) = y(n) + (h / 2) * (k1 + k2);
        x(n + 1) = x_next;

        if usa_exacta
            % paso local partiendo de la solucion exacta en x(n)
            y_ex_n = y_exact(x(n));
            y_ex_n1 = y_exact(x(n + 1));

            k1_ex = f(x(n), y_ex_n);
            y_pred_ex = y_ex_n + h * k1_ex;
            k2_ex = f(x(n + 1), y_pred_ex);

            y_loc_step = y_ex_n + (h / 2) * (k1_ex + k2_ex);

            e_local(n) = y_ex_n1 - y_loc_step;
            e_global(n + 1) = y_ex_n1 - y(n + 1);
        end
    end
end

