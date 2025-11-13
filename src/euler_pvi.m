function [x, y] = euler_pvi(f, x0, y0, h, N)
% euler_pvi  Resuelve un PVI de primer orden con el metodo de Euler explicito.
%
%   [x, y] = euler_pvi(f, x0, y0, h, N)
%
%   f  : handle a la funcion derivada f(x,y)
%   x0 : punto inicial
%   y0 : valor inicial y(x0)
%   h  : paso
%   N  : cantidad de pasos
%
%   x  : vector de puntos x_n
%   y  : vector de aproximaciones y_n

    x = zeros(N + 1, 1);
    y = zeros(N + 1, 1);

    x(1) = x0;
    y(1) = y0;

    for n = 1:N
        x(n + 1) = x(n) + h;
        y(n + 1) = y(n) + h * f(x(n), y(n));
    end
end

