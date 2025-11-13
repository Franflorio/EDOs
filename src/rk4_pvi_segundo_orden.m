function [x, y, yp] = rk4_pvi_segundo_orden(g, x0, y0, yp0, h, xfinal)
% rk4_pvi_segundo_orden  Resuelve un PVI de 2do orden
%   y'' = g(x, y, y') usando RK4 sobre el sistema equivalente.
%
%   [x, y, yp] = rk4_pvi_segundo_orden(g, x0, y0, yp0, h, xfinal)
%
%   g      : handle a g(x, y, yp) que define y'' = g(x, y, y')
%   x0     : punto inicial
%   y0     : valor inicial y(x0)
%   yp0    : valor inicial y'(x0)
%   h      : paso
%   xfinal : valor final de x
%
%   x      : vector de puntos x_n
%   y      : aproximaciones de y(x_n)
%   yp     : aproximaciones de y'(x_n)

    % definimos el sistema de 1er orden:
    % y1' = y2
    % y2' = g(x, y1, y2)
    f = @(x, Y) [Y(2);
                 g(x, Y(1), Y(2))];

    y0_vec = [y0; yp0];

    [x, Y] = rk4_sistema(f, x0, y0_vec, h, xfinal);

    y  = Y(:, 1);
    yp = Y(:, 2);
end

