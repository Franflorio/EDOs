function [x, Y] = euler_sistema(f, x0, y0, h, xfinal)
% euler_sistema  Metodo de Euler explicito para sistemas de EDO de primer orden.
%
%   [x, Y] = euler_sistema(f, x0, y0, h, xfinal)
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
%   Y      : matriz (N+1 x m), cada fila es el estado en x_n:
%            Y(n, :) = y_n'

    Nreal = (xfinal - x0) / h;
    N = round(Nreal);

    if abs(Nreal - N) > 1.0e-10
        warning('euler_sistema: El intervalo no es multiplo exacto de h. N = %d pasos.', N);
    end

    m = length(y0);

    x = zeros(N + 1, 1);
    Y = zeros(N + 1, m);

    x(1) = x0;
    Y(1, :) = y0(:)';

    for n = 1:N
        xn = x(n);
        yn = Y(n, :)';

        k1 = f(xn, yn);

        yn1 = yn + h * k1;

        x(n + 1) = xn + h;
        Y(n + 1, :) = yn1';
    end
end

