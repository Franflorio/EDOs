# Unidad 6 – Métodos Numéricos con EDO, Sistemas, 2º Orden y PVF

Repositorio de funciones en **MATLAB / Octave** para la Unidad 6 de Métodos Numéricos:

- Ecuaciones diferenciales ordinarias (EDO) de **primer orden** (PVI)
- **Sistemas** de EDO de primer orden
- EDO de **segundo orden** (vía sistemas)
- **Problemas de Valores de Frontera (PVF)**:
  - Método del **disparo**
  - Método de **diferencias finitas**

Además, el repositorio incluye un archivo `RESUMEN.md` con la **mini guía teórica** de la unidad.

---

## 1. PVI de primer orden (escalares)

EDO de la forma:

\[
y' = f(x,y), \quad y(x_0) = y_0
\]

### Funciones base

- **`euler_pvi.m`**  
  Método de **Euler explícito**.  
  ```matlab
  [x, y] = euler_pvi(f, x0, y0, h, xfinal)
  ```

- **`euler_mejorado_pvi.m`**  
  Método de **Euler mejorado / Heun** (orden 2).  
  ```matlab
  [x, y] = euler_mejorado_pvi(f, x0, y0, h, xfinal)
  ```

- **`rk2_pvi.m`**  
  **Runge–Kutta de orden 2** (punto medio).  
  ```matlab
  [x, y] = rk2_pvi(f, x0, y0, h, xfinal)
  ```

- **`rk4_pvi.m`**  
  **Runge–Kutta clásico de orden 4**.  
  ```matlab
  [x, y] = rk4_pvi(f, x0, y0, h, xfinal)
  ```

### Versiones con auditoría de error

- `euler_pvi_auditoria.m`  
- `euler_mejorado_pvi_auditoria.m`  
- `rk2_pvi_auditoria.m`  
- `rk4_pvi_auditoria.m`  

Firma genérica:

```matlab
[x, y, e_local, e_global] = metodo_pvi_auditoria(f, y_exact, x0, y0, h, xfinal)
```

- `y_exact` : handle a la solución exacta `y_exact(x)` (o `[]` si no se usa).
- `e_local` : error local de truncamiento aproximado en cada paso.
- `e_global`: error global `y_exact(x(i)) - y(i)` en cada nodo.

---

## 2. Sistemas de EDO de primer orden

Sistemas del tipo:

\[
\mathbf{y}' = \mathbf{F}(x, \mathbf{y}),\quad \mathbf{y}(x_0) = \mathbf{y}_0
\]

donde `y` es un vector (por ejemplo, `[y1; y2]`).

### Funciones

- **`euler_sistema.m`**  
  Euler explícito para sistemas.  
  ```matlab
  [x, Y] = euler_sistema(f, x0, y0, h, xfinal)
  ```

- **`euler_mejorado_sistema.m`**  
  Euler mejorado (Heun) para sistemas.  
  ```matlab
  [x, Y] = euler_mejorado_sistema(f, x0, y0, h, xfinal)
  ```

- **`rk4_sistema.m`**  
  Runge–Kutta clásico de orden 4 para sistemas.  
  ```matlab
  [x, Y] = rk4_sistema(f, x0, y0, h, xfinal)
  ```

Donde:

- `f(x, Y)` recibe un escalar `x` y un **vector columna** `Y`.
- `Y` es una matriz `(N+1) x m`:
  - `Y(:,1)` = aproximaciones de la primera variable,
  - `Y(:,2)` = de la segunda, etc.

Ejemplo de uso típico:

```matlab
% Sistema u' = (-2y - 6u)/120 ; y' = u
f = @(x, Y) [(-2 * Y(2) - 6 * Y(1)) / 120;
              Y(1)];

Y0 = [u0; y0];
[x, Y] = rk4_sistema(f, 0, Y0, 0.1, 10);
u = Y(:, 1);
y = Y(:, 2);
```

---

## 3. EDO de segundo orden (PVI)

Problemas de la forma:

\[
y'' = g(x,y,y'),\quad y(x_0)=y_0,\ y'(x_0)=y'_0
\]

Se resuelven convirtiendo a sistema de 1er orden:

\[
y_1 = y,\quad y_2 = y'
\Rightarrow
\begin{cases}
y_1' = y_2 \\
y_2' = g(x,y_1,y_2)
\end{cases}
\]

### Wrappers

- **`euler_pvi_segundo_orden.m`**  
- **`euler_mejorado_pvi_segundo_orden.m`**  
- **`rk4_pvi_segundo_orden.m`**  

Firma:

```matlab
[x, y, yp] = metodo_pvi_segundo_orden(g, x0, y0, yp0, h, xfinal)
```

- `g(x, y, yp)` define `y'' = g(x,y,y')`.
- `y`  : aproximaciones de `y(x)`.
- `yp` : aproximaciones de `y'(x)`.

---

## 4. PVF – Método del disparo

Problema de valores de frontera:

\[
y'' = g(x,y,y'),\quad a \le x \le b,
\]
\[
y(a) = \alpha,\quad y(b) = \beta
\]

Idea del método del disparo:

1. Introducir un parámetro \(s = y'(a)\) (pendiente inicial desconocida).
2. Para cada \(s\), resolver el PVI:
   - `y'' = g(x,y,y')`
   - `y(a) = alpha`, `y'(a) = s`
3. Definir la función disparo:
   \[
   \Phi(s) = y(b;s) - \beta
   \]
4. Usar un método de búsqueda de raíces (secante) para encontrar `s` tal que `Phi(s) ≈ 0`.

### Implementación

- **`disparo_rk4_segundo_orden.m`**  
  Usa RK4 de 2º orden (`rk4_pvi_segundo_orden`) + método de la secante sobre `s`.

Firma:

```matlab
[x, y, yp, s_opt, iter] = disparo_rk4_segundo_orden( ...
    g, a, b, alpha, beta, s1, s2, h, tol, maxit)
```

Parámetros:

- `g`     : función de la EDO `y'' = g(x,y,y')`
- `a, b`  : extremos del intervalo
- `alpha` : condición izquierda `y(a)`
- `beta`  : condición derecha `y(b)`
- `s1,s2` : dos disparos iniciales para `y'(a)`
- `h`     : paso de integración para RK4
- `tol`   : tolerancia para `|Phi(s)|`
- `maxit` : máximo de iteraciones de la secante

Devuelve:

- `x, y, yp` : solución final (aprox. de `y` y `y'`)
- `s_opt`    : pendiente inicial encontrada
- `iter`     : número de iteraciones usadas

---

## 5. PVF – Diferencias finitas (EDO lineales)

Para EDO lineales de 2º orden:

\[
y'' + p(x) y' + q(x) y = r(x),\quad a \le x \le b,
\]
\[
y(a) = \alpha,\quad y(b) = \beta
\]

Se discretiza el intervalo en una malla uniforme y se aproximan derivadas con diferencias finitas, obteniendo un **sistema lineal tridiagonal** en los nodos interiores.

### Implementaciones

- **`dif_finitas_pvf_lineal.m`**  
  Construye la malla, arma el sistema tridiagonal y calcula `y(x_i)`.

  ```matlab
  [x, y] = dif_finitas_pvf_lineal(p, q, rfun, a, b, alpha, beta, N)
  ```

  Donde:

  - `p(x), q(x), rfun(x)` son handles a las funciones del problema.
  - `N` es el número de subintervalos → `N+1` nodos.

- **`dif_finitas_pvf_lineal_auditoria.m`**  
  Versión con auditoría de error si se conoce la solución exacta `y_exact(x)`.

  ```matlab
  [x, y, e_global, tau] = dif_finitas_pvf_lineal_auditoria( ...
      p, q, rfun, a, b, alpha, beta, N, y_exact)
  ```

  Devuelve además:

  - `e_global(i) = y_exact(x(i)) - y(i)` → error global en cada nodo.
  - `tau(i)` → error de truncamiento (residual del esquema) en nodos interiores.

---

## 6. Resumen teórico: `RESUMEN.md`

En este mismo repositorio se incluye:

- **`RESUMEN.md`**  

Contiene una mini guía teórica de la Unidad 6:

- PVI de 1er orden y métodos de un paso (Euler, Euler mejorado, RK2, RK4)
- Sistemas de EDO y notación vectorial
- EDO de 2º orden y reducción a sistemas
- PVF:
  - Método del disparo
  - Método de diferencias finitas
- Sugerencias sobre **qué método usar cuándo**
- Checklist “modo parcial/TP” para PVI y PVF

Se recomienda leer primero `RESUMEN.md` para repasar teoría  
y luego explorar los `.m` según el tipo de problema a resolver.