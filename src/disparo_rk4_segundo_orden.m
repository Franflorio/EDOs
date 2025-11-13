function [x, y, yp, s_opt, iter] = disparo_rk4_segundo_orden(g, a, b, alpha, beta, s1, s2, h, tol, maxit)
% disparo_rk4_segundo_orden  Metodo del disparo para PVF de 2do orden
% usando RK4 y metodo de la secante sobre la pendiente inicial.
%
%   Problema:
%       y'' = g(x, y, y')    en el intervalo [a, b]
%       y(a) = alpha
%       y(b) = beta
%
%   Idea del metodo del disparo:
%       - No conocemos y'(a), la llamamos s (pendiente inicial).
%       - Para cada s, resolvemos el PVI:
%             y'' = g(x, y, y'),  y(a) = alpha,  y'(a) = s
%         y obtenemos un valor aproximado de y(b; s).
%       - Definimos Phi(s) = y(b; s) - beta.
%       - Buscamos s tal que Phi(s) = 0 (la solucion cumple la condicion en b).
%       - Aqui usamos el metodo de la secante para encontrar esa s.
%
%   [x, y, yp, s_opt, iter] =
%       disparo_rk4_segundo_orden(g, a, b, alpha, beta, s1, s2, h, tol, maxit)
%
%   Entradas:
%       g     : handle a g(x, y, yp) tal que y'' = g(x, y, y')
%       a, b  : extremos del intervalo [a, b]
%       alpha : condicion de frontera en a: y(a) = alpha
%       beta  : condicion de frontera en b: y(b) = beta
%       s1    : primer valor de prueba para la pendiente inicial y'(a)
%       s2    : segundo valor de prueba para la pendiente inicial y'(a)
%       h     : paso de integracion para RK4 (en x)
%       tol   : tolerancia para |Phi(s)|
%       maxit : maximo numero de iteraciones de la secante
%
%   Salidas:
%       x     : vector de puntos x_n en [a, b]
%       y     : aproximaciones de y(x_n) para la pendiente s_opt
%       yp    : aproximaciones de y'(x_n) para la pendiente s_opt
%       s_opt : pendiente inicial y'(a) encontrada por el metodo del disparo
%       iter  : cantidad de iteraciones realizadas por la secante
%
%   Nota:
%       Esta funcion llama internamente a:
%           rk4_pvi_segundo_orden(g, a, alpha, s, h, b)
%       que debe estar definida por separado (usa RK4 para el PVI de 2do orden).

    % Evaluamos Phi(s1) y Phi(s2).
    % Cada Phi(si) es la diferencia entre y(b; si) y beta.
    Phi1 = phi_disparo(g, a, b, alpha, beta, s1, h);
    Phi2 = phi_disparo(g, a, b, alpha, beta, s2, h);

    % Si Phi1 y Phi2 tienen el mismo signo, en principio no garantizamos
    % que haya una raiz en el intervalo [s1, s2]. Igual avisamos, pero
    % el metodo de la secante puede llegar a converger de todas formas.
    if Phi1 * Phi2 > 0
        warning('disparo_rk4_segundo_orden: Phi(s1) y Phi(s2) tienen el mismo signo.');
    end

    iter = 0;

    % Bucle principal del metodo de la secante para encontrar s:
    %   s_{k+1} = s_k - Phi(s_k) * (s_k - s_{k-1}) / (Phi(s_k) - Phi(s_{k-1}))
    for k = 1:maxit
        iter = k;

        % Si Phi2 - Phi1 es muy pequeño, el metodo de la secante puede
        % explotar numericamente. Cortamos la iteracion.
        if abs(Phi2 - Phi1) < 1.0e-14
            warning('disparo_rk4_segundo_orden: Diferencia Phi2 - Phi1 casi nula. Se detiene.');
            break;
        end

        % Paso de la secante:
        %   s_new es la nueva aproximacion de la pendiente inicial.
        s_new = s2 - Phi2 * (s2 - s1) / (Phi2 - Phi1);

        % Calculamos Phi(s_new) integrando de nuevo el PVI
        % con esa pendiente inicial.
        Phi_new = phi_disparo(g, a, b, alpha, beta, s_new, h);

        % Si ya estamos suficientemente cerca de la condicion en b,
        % (es decir, |Phi(s_new)| < tol), aceptamos s_new y salimos.
        if abs(Phi_new) < tol
            s2 = s_new;
            Phi2 = Phi_new;
            break;
        end

        % Actualizamos los valores para la siguiente iteracion:
        % El viejo "s2" pasa a ser "s1" y s_new pasa a ser el nuevo "s2".
        s1 = s2;
        Phi1 = Phi2;

        s2 = s_new;
        Phi2 = Phi_new;
    end

    % s_opt es la pendiente inicial final encontrada por la secante.
    s_opt = s2;

    % Con la pendiente inicial s_opt ya fijada, integramos una ultima vez
    % el PVI en todo [a, b] para obtener las soluciones finales y(x), y'(x).
    [x, y, yp] = rk4_pvi_segundo_orden(g, a, alpha, s_opt, h, b);
end

function Phi = phi_disparo(g, a, b, alpha, beta, s, h)
% phi_disparo  Calcula la funcion Phi(s) = y(b; s) - beta
% para el metodo del disparo.
%
%   y'' = g(x, y, y')
%   y(a) = alpha
%   y'(a) = s
%
%   La idea:
%       - Resolvemoss el PVI en [a, b] con condiciones iniciales
%         y(a) = alpha, y'(a) = s.
%       - Tomamos el valor aproximado y(b; s) (ultimo punto del vector y).
%       - Phi(s) = y(b; s) - beta indica cuan lejos estamos de cumplir
%         la condicion de frontera y(b) = beta.

    % Resolvemos el PVI de 2do orden con RK4 para este valor de s.
    [x, y, yp] = rk4_pvi_segundo_orden(g, a, alpha, s, h, b);

    % El ultimo valor de y corresponde a x = b (o muy cercano).
    yb = y(end);

    % Phi(s) mide el error en la condicion de frontera en b.
    Phi = yb - beta;
end

