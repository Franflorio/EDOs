function [x, y, e_global, tau] = dif_finitas_pvf_lineal_auditoria(p, q, rfun, a, b, alpha, beta, N, y_exact)
% dif_finitas_pvf_lineal_auditoria
%   Metodo de diferencias finitas para PVF lineal de 2do orden con auditoria.
%
%   EDO:
%       y'' + p(x) y' + q(x) y = r(x),   a <= x <= b
%       y(a) = alpha
%       y(b) = beta
%
%   [x, y, e_global, tau] =
%       dif_finitas_pvf_lineal_auditoria(p, q, rfun, a, b, alpha, beta, N, y_exact)
%
%   Entradas:
%       p       : handle a p(x)
%       q       : handle a q(x)
%       rfun    : handle a r(x)
%       a, b    : extremos del intervalo
%       alpha   : condicion de frontera en a, y(a) = alpha
%       beta    : condicion de frontera en b, y(b) = beta
%       N       : numero de subintervalos (habra N+1 nodos)
%       y_exact : handle a la solucion exacta y(x).
%                 Si es [] o no es function_handle, no se calculan errores
%
%   Salidas:
%       x        : vector columna de nodos x_i, i = 0..N
%       y        : aproximaciones de y(x_i) en cada nodo
%       e_global : error global en cada nodo, y_exact(x_i) - y_i.
%                  Vector vacio si no se pasa y_exact.
%       tau      : error de truncamiento (local) en cada nodo interior.
%                  tau(i) es la residual del esquema en x_i usando
%                  la solucion exacta. En los extremos se pone NaN.
%                  Vector vacio si no se pasa y_exact.

    % Paso de malla
    h = (b - a) / N;

    % Nodos x_0,...,x_N
    x = linspace(a, b, N + 1)';

    % Vector solucion (incluye fronteras)
    y = zeros(N + 1, 1);
    y(1) = alpha;       % y_0 = alpha
    y(end) = beta;      % y_N = beta

    % Cantidad de incognitas interiores: y_1,...,y_{N-1}
    m = N - 1;

    % Matriz tridiagonal y vector de terminos independientes
    A = zeros(m, m);
    d = zeros(m, 1);

    % Construccion del sistema para los nodos interiores
    for j = 1:m
        % Nodo global i (entre 1 y N-1, pero recorda que en MATLAB
        % x(1) es x_0, entonces x(i) es x_{i-1} si usaramos otro indice;
        % aqui definimos i = j+1 para que x(i) sea x_j global).
        i = j + 1;          % i recorre 2..N
        xi = x(i);

        pi = p(xi);
        qi = q(xi);
        ri = rfun(xi);

        % Coeficientes del esquema:
        % A_i y_{i-1} + B_i y_i + C_i y_{i+1} = r_i h^2
        Ai = 1 - (pi * h) / 2;
        Bi = -2 + qi * h * h;
        Ci = 1 + (pi * h) / 2;
        di = ri * h * h;

        % Tratamiento de la frontera izquierda: y_0 = alpha
        if j == 1
            % Aparece A_1 * y_0, que es conocido, se pasa al lado derecho.
            di = di - Ai * alpha;
        else
            % A_i multiplica a y_{i-1}, que es una incognita interior
            A(j, j - 1) = Ai;
        end

        % Tratamiento de la frontera derecha: y_N = beta
        if j == m
            % Aparece C_{N-1} * y_N, conocido, se pasa al lado derecho.
            di = di - Ci * beta;
        else
            % C_i multiplica a y_{i+1}, que es incognita interior
            A(j, j + 1) = Ci;
        end

        % Coeficiente diagonal B_i para y_i interior
        A(j, j) = Bi;
        d(j) = di;
    end

    % Resolver el sistema lineal para las incognitas interiores
    y_int = A \ d;

    % Insertar en el vector completo de solucion
    y(2:N) = y_int;

    % Auditoria: errores solo si se paso una solucion exacta
    if isa(y_exact, 'function_handle')
        % Error global en cada nodo: y_exact(x_i) - y_i
        y_ex = y_exact(x);
        e_global = y_ex - y;

        % Error de truncamiento (local) en nodos interiores:
        % tau(i) = A_i y_ex(i-1) + B_i y_ex(i) + C_i y_ex(i+1) - r_i h^2
        tau = NaN(size(y));   % extremos quedan como NaN

        for j = 1:m
            i = j + 1;
            xi = x(i);

            pi = p(xi);
            qi = q(xi);
            ri = rfun(xi);

            Ai = 1 - (pi * h) / 2;
            Bi = -2 + qi * h * h;
            Ci = 1 + (pi * h) / 2;

            yi_minus = y_ex(i - 1);
            yi       = y_ex(i);
            yi_plus  = y_ex(i + 1);

            tau(i) = Ai * yi_minus + Bi * yi + Ci * yi_plus - ri * h * h;
        end
    else
        e_global = [];
        tau = [];
    end
end

