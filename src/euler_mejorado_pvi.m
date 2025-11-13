function [x, y] = euler_mejorado_pvi(f, x0, y0, h, xfinal)
% euler_mejorado_pvi  Resuelve un PVI de primer orden con el metodo de Euler mejorado (Heun).
%
%   [x, y] = euler_mejorado_pvi(f, x0, y0, h, xfinal)
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
        warning('El intervalo no es multiplo exacto de h. Se usa N = %d pasos.', N);
    end

    x = zeros(N + 1, 1);
    y = zeros(N + 1, 1);

    x(1) = x0;
    y(1) = y0;

    for n = 1:N
        % pendiente en el inicio
        k1 = f(x(n), y(n));

        % prediccion con Euler simple
        y_pred = y(n) + h * k1;
        x_next = x(n) + h;

        % pendiente en el final usando la prediccion
        k2 = f(x_next, y_pred);

        % correccion con promedio de pendientes
        y(n + 1) = y(n) + (h / 2) * (k1 + k2);
        x(n + 1) = x_next;
    end
end

