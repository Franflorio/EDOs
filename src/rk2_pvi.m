function [x, y] = rk2_pvi(f, x0, y0, h, xfinal)
% rk2_pvi  Resuelve un PVI de primer orden con Runge Kutta de orden 2
%          (metodo del punto medio).
%
%   [x, y] = rk2_pvi(f, x0, y0, h, xfinal)
%
%   f      : handle a la funcion derivada f(x, y)
%   x0     : punto inicial
%   y0     : valor inicial y(x0)
%   h      : paso
%   xfinal : valor final de x
%
%   x      : vector de puntos x_n
%   y      : vector de aproximaciones y_n

    % calcular cantidad de pasos
    Nreal = (xfinal - x0) / h;
    N = round(Nreal);

    if abs(Nreal - N) > 1.0e-10
        warning('rk2_pvi: El intervalo no es multiplo exacto de h. N = %d pasos.', N);
    end

    x = zeros(N + 1, 1);
    y = zeros(N + 1, 1);

    x(1) = x0;
    y(1) = y0;

    for n = 1:N
        k1 = f(x(n), y(n));
        x_med = x(n) + h * 0.5;
        y_med = y(n) + h * 0.5 * k1;

        k2 = f(x_med, y_med);

        x(n + 1) = x(n) + h;
        y(n + 1) = y(n) + h * k2;
    end
end

