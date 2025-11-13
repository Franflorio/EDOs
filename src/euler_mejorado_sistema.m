function [x, Y] = euler_mejorado_sistema(f, x0, y0, h, xfinal)
% euler_mejorado_sistema  Metodo de Euler mejorado (Heun)
%                         para sistemas de EDO de primer orden.
%
%   [x, Y] = euler_mejorado_sistema(f, x0, y0, h, xfinal)
%
%   f      : handle a la funcion derivada f(x, y), y vector columna
%   x0     : punto inicial
%   y0     : vector columna con las condiciones iniciales
%   h      : paso
%   xfinal : valor final de x
%
%   x      : vector de puntos x_n
%   Y      : matriz (N+1 x m), cada fila es y_n'

    Nreal = (xfinal - x0) / h;
    N = round(Nreal);

    if abs(Nreal - N) > 1.0e-10
        warning('euler_mejorado_sistema: El intervalo no es multiplo exacto de h. N = %d pasos.', N);
    end

    m = length(y0);

    x = zeros(N + 1, 1);
    Y = zeros(N + 1, m);

    x(1) = x0;
    Y(1, :) = y0(:)';

    for n = 1:N
        xn = x(n);
        yn = Y(n, :)';

        % pendiente en el inicio
        k1 = f(xn, yn);

        % prediccion con Euler simple
        x_next = xn + h;
        y_pred = yn + h * k1;

        % pendiente en el final
        k2 = f(x_next, y_pred);

        % correccion con promedio de pendientes
        yn1 = yn + (h / 2) * (k1 + k2);

        x(n + 1) = x_next;
        Y(n + 1, :) = yn1';
    end
end

