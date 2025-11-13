 % euler_mejorado_clase.m
 %Script sacado del pdf ejemplo de ecuaciones lineales.
 M=276.35;
 K=28000.0;
 c=1750.0;
  f1 = @(x, y, u) (28000 * y - 1750.0 * u) / 276.35;
  f2 = @(x, y, u) u;
 h = 0.05;
 x = 0;
 y = 9;
 u = 0;
 iteraciones = round(0.6/h);
 for i = 1:iteraciones
 x(i+1) = x(i) + h;
 y_euler = y(i) + h*f2(x(i),y(i),u(i));
 u_euler = u(i) + h*f1(x(i),y(i),u(i));
 y(i+1) = y(i) + h/2 * (f2(x(i),y(i),u(i)) + f2(x(i+1), y_euler, u_euler));
 u(i+1) = u(i) + h/2 * (f1(x(i),y(i),u(i)) + f1(x(i+1), y_euler, u_euler));
 end
