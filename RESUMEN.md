# Unidad 6 – EDO, PVI y PVF (Mini guía)

## 0. Panorama de la unidad

### Tipos de problemas

1. **PVI de 1er orden**  
   \[
   y' = f(x,y),\quad y(x_0) = y_0
   \]

2. **Sistemas de EDO de 1er orden**  
   \[
   \mathbf{y}' = \mathbf{F}(x,\mathbf{y}),\quad \mathbf{y}(x_0) = \mathbf{y}_0
   \]

3. **EDO de 2º orden (PVI)**  
   \[
   y'' = g(x,y,y'),\quad y(x_0)=y_0,\ y'(x_0)=y'_0
   \]

4. **PVF (Problemas de Valores de Frontera)**  
   \[
   y'' = g(x,y,y'),\ a\le x\le b,\quad y(a)=\alpha,\ y(b)=\beta
   \]

---

## 1. PVI de primer orden

### Modelo

- EDO:  
  \[
  y' = f(x,y)
  \]
- Dato: \(y(x_0)=y_0\)
- Objetivo: aproximar \(y(x)\) en una malla \(x_0, x_1,\dots,x_N\) con paso \(h\).

### Métodos de un paso (idea general)

Todos tienen la forma:

\[
y_{n+1} = y_n + h\;\Phi(f, x_n, y_n, h)
\]

### Métodos vistos

| Método            | Fórmula base                                               | Orden global | Evaluaciones de \(f\) por paso |
|------------------|------------------------------------------------------------|--------------|--------------------------------|
| Euler            | \(y_{n+1} = y_n + h f(x_n,y_n)\)                           | 1            | 1                              |
| Euler mejorado   | promedio de pendientes (inicio y final)                    | 2            | 2                              |
| RK2 (punto medio)| usa la pendiente en el punto medio                         | 2            | 2                              |
| RK4 clásico      | combinación \(k_1, k_2, k_3, k_4\)                         | 4            | 4                              |

- **Euler**: simple, barato, pero poco preciso (orden 1).  
- **Euler mejorado / RK2**: mejor precisión con poco esfuerzo extra.  
- **RK4**: estándar para buena precisión en muchos problemas.

---

## 2. Sistemas de EDO de primer orden

### Modelo

Ejemplo con dos ecuaciones:

\[
\begin{cases}
y_1' = f_1(x,y_1,y_2) \\
y_2' = f_2(x,y_1,y_2)
\end{cases}
\]

Se agrupan en un vector:

\[
\mathbf{y} =
\begin{pmatrix}
y_1 \\ y_2
\end{pmatrix},
\quad
\mathbf{y}' = \mathbf{F}(x,\mathbf{y})
\]

### Idea clave

> Los mismos métodos (Euler, RK2, RK4) se aplican igual,  
> solo que ahora `y` y los `k` son **vectores**.

En cada paso, por ejemplo en RK4:

\[
\begin{aligned}
\mathbf{k}_1 &= \mathbf{F}(x_n, \mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{F}\!\left(x_n + \tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{F}\!\left(x_n + \tfrac{h}{2},\; \mathbf{y}_n + \tfrac{h}{2}\mathbf{k}_2\right) \\
\mathbf{k}_4 &= \mathbf{F}\!\left(x_n + h,\; \mathbf{y}_n + h\mathbf{k}_3\right) \\
\mathbf{y}_{n+1} &= \mathbf{y}_n + \tfrac{h}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)
\end{aligned}
\]

Cada componente de \(\mathbf{y}\) representa una variable del sistema (posición, velocidad, corrientes, etc.).

---

## 3. EDO de 2º orden (PVI) → sistema de 1er orden

### Modelo

\[
y'' = g(x,y,y'),\quad y(x_0)=y_0,\ y'(x_0)=y'_0
\]

### Truco estándar de reducción de orden

Definimos:

\[
y_1 = y,\quad y_2 = y'
\]

Entonces:

\[
\begin{cases}
y_1' = y_2 \\
y_2' = g(x,y_1,y_2)
\end{cases}
\]

Esto es un sistema de 2 EDO de 1er orden.  
Se puede resolver con cualquier método para sistemas (Euler, RK2, RK4).

### Ventaja

- Toda la teoría y el código para sistemas se reutiliza.
- Para el usuario, es cómodo tener funciones “wrapper” que manejen la conversión y devuelvan directamente `y` y `y'`.

---

## 4. PVF: PVI vs PVF y métodos

### Diferencia PVI vs PVF

- **PVI** (problema de valor inicial):  
  todos los datos están en un punto:
  \[
  y(x_0)=y_0, \quad y'(x_0)=y'_0
  \]

- **PVF** (problema de valores de frontera):  
  las condiciones están en los extremos:
  \[
  y(a) = \alpha,\quad y(b) = \beta
  \]

En PVF no conocés \(y'(a)\) → no podés aplicar directamente los métodos de PVI.

Dos métodos principales:

1. **Método del disparo**  
2. **Diferencias finitas**

---

## 4.1. Método del disparo

### Idea conceptual

1. Introducís un parámetro \(s = y'(a)\) (pendiente inicial desconocida).
2. Para cada \(s\) resolvés el PVI:
   \[
   y'' = g(x,y,y'),\quad y(a)=\alpha,\ y'(a)=s
   \]
   obteniendo una solución \(y(x;s)\).
3. En el extremo derecho calculás:
   \[
   \Phi(s) = y(b;s) - \beta
   \]
4. Querés que se cumpla \(\Phi(s)=0\), es decir, que la solución numérica respete la condición de frontera en \(b\).

### Resolución numérica

Buscar la raíz de \(\Phi(s)\):

- **Bisección**: si se encuentra un intervalo \([s_1,s_2]\) con \(\Phi(s_1)\Phi(s_2)<0\).
- **Secante**: actualización:
  \[
  s_{\text{nuevo}} = s_2 - \Phi(s_2)\,\frac{s_2 - s_1}{\Phi(s_2) - \Phi(s_1)}
  \]

En cada evaluación de \(\Phi(s)\) hay que resolver **un PVI completo**.

### Ventajas

- Reutiliza por completo los métodos para PVI (RK4, etc.).
- Intuitivo: “disparar” con diferentes pendientes hasta “pegarle” a la frontera derecha.

### Desventajas

- Puede ser inestable o complicado si el problema es muy sensible a las condiciones iniciales.
- Requiere resolver el PVI varias veces.

---

## 4.2. Diferencias finitas para PVF lineales

### Modelo lineal

\[
y'' + p(x) y' + q(x) y = r(x),\quad a\le x\le b,
\]
\[
y(a)=\alpha,\ y(b)=\beta.
\]

### Pasos del método

1. **Malla uniforme**:  
   \[
   x_i = a + i h,\quad h = \frac{b-a}{N},\ i=0,\dots,N
   \]

2. **Aproximaciones en diferencias**:
   \[
   y''(x_i) \approx \frac{y_{i+1} - 2y_i + y_{i-1}}{h^2},\quad
   y'(x_i) \approx \frac{y_{i+1} - y_{i-1}}{2h}
   \]

3. **Ecuación en cada nodo interior** \(x_i\), \(i=1,\dots,N-1\):

   \[
   A_i y_{i-1} + B_i y_i + C_i y_{i+1} = r_i h^2,
   \]
   donde:
   \[
   \begin{aligned}
   A_i &= 1 - \frac{p_i h}{2} \\
   B_i &= -2 + q_i h^2 \\
   C_i &= 1 + \frac{p_i h}{2} \\
   p_i &= p(x_i),\ q_i = q(x_i),\ r_i = r(x_i)
   \end{aligned}
   \]

4. **Incorporar condiciones de frontera**:
   - \(y_0 = \alpha\)
   - \(y_N = \beta\)

   Su efecto es mover términos con \(y_0,y_N\) al lado derecho del sistema.

5. **Sistema tridiagonal** en las incógnitas \(y_1,\dots,y_{N-1}\).  
   Se resuelve con álgebra lineal (backslash, métodos específicos, etc.).

### Ventajas

- Muy natural para EDO **lineales**.
- Lleva a un sistema lineal que se puede resolver de una vez.
- Matriz tridiagonal → buena eficiencia numérica.

### Desventajas

- Hay que definir una malla y tamaño de paso `h` adecuados.
- Para EDO no lineales, el sistema resultante deja de ser lineal.

---

## 5. ¿Qué método usar cuándo?

### PVI (1er o 2º orden)

- **Euler**:
  - Uso pedagógico o problemas muy simples.
  - Sirve para entender la dinámica, no es muy preciso.

- **Euler mejorado / RK2**:
  - Buena opción cuando querés mejor precisión que Euler sin mucho costo extra.

- **RK4**:
  - Opción estándar para buena precisión en muchos problemas.
  - Buen balance entre costo y exactitud.

- **ode45** (MATLAB/Octave):
  - Para trabajo práctico “real”: método adaptativo RK(4,5).
  - Ajusta el paso automáticamente según tolerancias de error.

### PVF

- **Disparo**:
  - Buena opción cuando ya tenés implementado un buen solver de PVI (RK4).
  - Adecuado para EDO potencialmente no lineales.
  - Intuitivo pero requiere varias integraciones del PVI.

- **Diferencias finitas**:
  - Ideal para EDO **lineales** de 2º orden.
  - Se formula una sola vez como sistema lineal.
  - Muy común en parciales: “discretizar y plantear el sistema”.

---

## 6. Checklist rápido (modo parcial / TP)

### Si el problema es un PVI

1. Identificar:
   - ¿Es de 1er orden? ¿de 2º?  
2. Si es de 2º orden:
   - Reducir a sistema de 1er orden o usar wrappers.
3. Elegir método:
   - Euler, Euler mejorado, RK2, RK4 (según te pidan).
4. Armar malla:
   - `x0, x1, ..., xN` con paso `h`.
5. Aplicar el esquema paso a paso y construir tabla de valores.
6. Si existe solución exacta:
   - Calcular error global en cada nodo.
   - Comentar comportamiento del error al variar `h`.

### Si el problema es un PVF

1. Ver qué método te piden:
   - **Disparo**:
     - Plantear PVI con `y'(a)=s`.
     - Definir \(\Phi(s) = y(b;s) - \beta\).
     - Mostrar al menos dos disparos (`s1`, `s2`) y una interpolación tipo secante.
   - **Diferencias finitas**:
     - Escribir malla y fórmulas en diferencias.
     - Derivar ecuaciones para cada nodo interior.
     - Mostrar coeficientes \(A_i,B_i,C_i\).
     - Indicar cómo se incorporan las condiciones de frontera.
     - Plantear el sistema lineal resultante (y resolver si te lo piden).
