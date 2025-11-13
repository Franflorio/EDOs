function [x, y] = dif_finitas_pvf_lineal(p, q, rfun, a, b, alpha, beta, N)
% dif_finitas_pvf_lineal  Metodo de diferencias finitas para PVF lineal
% de 2do orden:
%
%   y'' + p(x) y' + q(x) y = r(x),   a <= x <= b
%   y(a) = alpha
%   y(b) = beta
%
%   [x, y] = dif_finitas_pvf_lineal(p, q, rfun, a, b, alpha, beta, N)
%
%   Entradas:
%       p     : handle a p(x)
%       q     : handle a q(x)
%       rfun  : handle a r(x)
%       a, b  : extremos del intervalo
%       alpha : condicion de frontera en a, y(a) = alpha
%       beta  : condicion de frontera en b, y(b) = beta
%       N     : numero de subintervalos (habra N+1 nodos)
%
%   Salidas:
%       x     : vector columna de nodos x_i, i = 0..N
%       y     : aproximaciones de y(x_i) en cada nodo
%
%   Nota:
%       - Se utiliza malla uniforme con paso h = (b - a) / N.
%       - Se aproxima:
%           y''(x_i) ~ (y_{i+1} - 2 y_i + y_{i-1}) / h^2
%           y'(x_i)  ~ (y_{i+1} - y_{i-1}) / (2 h)
%       - El sistema resultante es tridiagonal en las incognitas
%         y_1, ..., y_{N-1}.

    % Paso de malla
    h = (b - a) / N;

    % Nodos x_0,...,x_N
    x = linspace(a, b, N + 1)';

    % Inicializamos vector de solucion completo, incluyendo fronteras
    y = zeros(N + 1, 1);
    y(1) = alpha;       % y_0 = alpha
    y(end) = beta;      % y_N = beta

    % Cantidad de puntos interiores (incognitas)
    m = N - 1;          % y_1,...,y_{N-1}

    % Matriz tridiagonal y vector de terminos independientes
    A = zeros(m, m);
    d = zeros(m, 1);

    % Recorremos nodos interiores: i = 1..N-1
    % j = i en "indice interior" (1..m) corresponde al nodo x_i
    for j = 1:m
        i = j + 1;              % nodo global x_i (2..N)
        xi = x(i);

        pi = p(xi);
        qi = q(xi);
        ri = rfun(xi);

        % Coeficientes A_i, B_i, C_i de la ecuacion:
        % A_i y_{i-1} + B_i y_i + C_i y_{i+1} = r_i h^2
        Ai = 1 - (pi * h) / 2;
        Bi = -2 + qi * h * h;
        Ci = 1 + (pi * h) / 2;
        di = ri * h * h;

        % Tratamiento de A_i y Ci segun sea el primer o ultimo interior
        % Recordar:
        %   y_0 = alpha, y_N = beta son conocidos.
        % Para i=1 (j=1) aparece y_0 -> pasa al lado derecho.
        if j == 1
            % A_1 * y_0 se pasa restando al lado derecho
            di = di - Ai * alpha;
        else
            % A_i multiplica a y_{i-1}, que es una incognita interior
            A(j, j - 1) = Ai;
        end

        % Para i=N-1 (j=m) aparece y_N -> pasa al lado derecho.
        if j == m
            % C_{N-1} * y_N se pasa restando
            di = di - Ci * beta;
        else
            % C_i multiplica a y_{i+1}, incognita interior
            A(j, j + 1) = Ci;
        end

        % Termino diagonal siempre corresponde a y_i interior
        A(j, j) = Bi;
        d(j) = di;
    end

    % Resolvemos el sistema tridiagonal A * y_int = d
    y_int = A \ d;

    % Insertamos las soluciones interiores en el vector completo
    y(2:N) = y_int;
end

