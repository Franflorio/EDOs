function [t, Y] = euler_vectorizado(f, tspan, y0, h)
    % f     : función que representa el sistema de EDOs, f(t, y)
    % tspan : intervalo de tiempo [t0, tf]
    % y0    : vector columna con las condiciones iniciales
    % h     : tamaño de paso

    t0 = tspan(1);    % tiempo inicial
    tf = tspan(2);    % tiempo final
    t = t0:h:tf;      % vector de tiempos

    N = length(t);    % número de pasos de tiempo
    m = length(y0);   % tamaño del sistema de ecuaciones

    Y = zeros(m, N);  % matriz para almacenar las soluciones
    Y(:, 1) = y0;     % condiciones iniciales

    % Método de Euler vectorizado
    for n = 1:N-1
        Y(:, n+1) = Y(:, n) + h * f(t(n), Y(:, n));
    end
end
