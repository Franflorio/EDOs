function [x, y, yp] = euler_mejorado_pvi_segundo_orden(g, x0, y0, yp0, h, xfinal)
% euler_mejorado_pvi_segundo_orden  Euler mejorado (Heun) para
%   PVI de 2do orden y'' = g(x, y, y').
%
%   [x, y, yp] = euler_mejorado_pvi_segundo_orden(g, x0, y0, yp0, h, xfinal)
%
%   g      : handle a g(x, y, yp)
%   x0     : punto inicial
%   y0     : valor inicial y(x0)
%   yp0    : valor inicial y'(x0)
%   h      : paso
%   xfinal : valor final de x
%
%   x      : vector de puntos x_n
%   y      : aproximaciones de y(x_n)
%   yp     : aproximaciones de y'(x_n)

    f = @(x, Y) [Y(2);
                 g(x, Y(1), Y(2))];

    y0_vec = [y0; yp0];

    [x, Y] = euler_mejorado_sistema(f, x0, y0_vec, h, xfinal);

    y  = Y(:, 1);
    yp = Y(:, 2);
end

