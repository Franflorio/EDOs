function [x, Y] = rk4_sistema(f, x0, y0, h, xfinal)
% rk4_sistema  Runge Kutta clasico de orden 4 para sistemas de EDO de primer orden.
%
%   [x, Y] = rk4_sistema(f, x0, y0, h, xfinal)
%
%   f      : handle a la funcion derivada f(x, y), donde
%            y es un vector columna y f devuelve un vector columna
%            del mismo tamaño.
%   x0     : punto inicial
%   y0     : vector columna con las condiciones iniciales
%   h      : paso
%   xfinal : valor final de x
%
%   x      : vector de puntos x_n (N+1 x 1)
%   Y      : matriz (N+1 x m), donde cada fila es y_n' (estado en x_n)

    Nreal = (xfinal - x0) / h;
    N = round(Nreal);

    if abs(Nreal - N) > 1.0e-10
        warning('rk4_sistema: El intervalo no es multiplo exacto de h. N = %d pasos.', N);
    end

    m = length(y0);          % dimension del sistema

    x = zeros(N + 1, 1);
    Y = zeros(N + 1, m);

    x(1) = x0;
    Y(1, :) = y0(:)';        % guardamos como fila

    for n = 1:N
        yn = Y(n, :)';
        xn = x(n);

        k1 = f(xn, yn);
        k2 = f(xn + h * 0.5, yn + h * 0.5 * k1);
        k3 = f(xn + h * 0.5, yn + h * 0.5 * k2);
        k4 = f(xn + h,       yn + h * k3);

        yn1 = yn + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

        x(n + 1) = xn + h;
        Y(n + 1, :) = yn1';
    end
end

