function [x, y, e_local, e_global] = euler_pvi_auditoria(f, y_exact, x0, y0, h, xfinal)
% euler_pvi_auditoria  Metodo de Euler con auditoria de errores.
%
%   [x, y, e_local, e_global] = euler_pvi_auditoria(f, y_exact, x0, y0, h, xfinal)
%
%   f       : handle a la funcion derivada f(x, y)
%   y_exact : handle a la solucion exacta y(x). Si es [] no se calculan errores.
%   x0, y0  : condicion inicial y(x0) = y0
%   h       : paso
%   xfinal  : valor final de x
%
%   x, y       : aproximaciones numericas
%   e_local    : error local por paso (y_exact(x_{n+1}) - y_loc_step)
%   e_global   : error global en cada punto (y_exact(x_n) - y_n)

    Nreal = (xfinal - x0) / h;
    N = round(Nreal);

    if abs(Nreal - N) > 1.0e-10
        warning('euler_pvi_auditoria: El intervalo no es multiplo exacto de h. N = %d pasos.', N);
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
        x(n + 1) = x(n) + h;
        y(n + 1) = y(n) + h * f(x(n), y(n));

        if usa_exacta
            y_ex_n = y_exact(x(n));
            y_ex_n1 = y_exact(x(n + 1));

            % paso local de Euler partiendo de la solucion exacta
            y_loc_step = y_ex_n + h * f(x(n), y_ex_n);

            e_local(n) = y_ex_n1 - y_loc_step;
            e_global(n + 1) = y_ex_n1 - y(n + 1);
        end
    end
end

