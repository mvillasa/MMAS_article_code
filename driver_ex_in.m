% author: Saul Buitrago (USB)
% script manejador del optimizador
% ejemplos:
%---- funcion 1
%eval_f = 10.*(z(2)-z(1)*z(1))^2 + (1.-z(1))^2  en [-1.2,0.8]x[-1.2,0.8]
%---- funcion 2
%eval_f = (z(1)-1.)^2+(z(2)-1.)^2-cos(18*(z(2)-1.))-cos(18*(z(1)-1.)) en [0,2]x[0,2]
%---- funcion 3
%eval_funcion = abs(-z(1)+z(2)-2.)+2. en [-2,0]x[0,2]
% datos
poblacion=20; input_no=2; cycle = 80; example = 1;
% grafico de la funcion
switch example
  case 3
    x1=-2.; x2=0.; y1=0.; y2=2.; % funcion 3  [-2,0]x[0,2]
    x = linspace(x1,x2,50); y = linspace(y1,y2,50);
  case 2
    x1=0.; x2=2.; y1=0.; y2=2.; % funcion 2  [0,2]x[0,2]
    x = linspace(x1,x2,100); y = x;
  case 1
    x1=-1.2; x2=0.8; y1=-1.2; y2=0.8; % funcion 1  [-1.2,0.8]x[-1.2,0.8]
    x = linspace(x1,x2,100); y = linspace(y1,y2,100);
end
figure(1);
%x = linspace(x1,x2,100); y = x;
[xx,yy] = meshgrid(x,y);
switch example
  case 3
    zz = abs(-xx+yy-2)+2; % funcion 3  [-2,0]x[0,2]
  case 2
    zz = (xx-1).^2+(yy-1).^2-cos(18*(xx-1))-cos(18*(yy-1)); % funcion 2  [0,2]x[0,2]
  case 1
    zz = (10*(yy-xx.^2)).^2+(1-xx).^2; % funcion 1  [-1.2,0.8]x[-1.2,0.8]
end
subplot(1, 2, 1), mesh(xx,yy,zz);
subplot(1, 2, 2), contour(x,y,zz,20);
% llamado al optimizador
[xx,fxy,deltf10,varianza,icont5] = ...
       opt_ex_in(poblacion,input_no,cycle,example);
disp('****Resultados en el main****');
parametros = [xx(1),xx(2)];
fprintf('funcion objetivo = %f\n',fxy);
disp('parametros calculados');
fprintf('x1= %f - x2= %f\n',...
    parametros(1),parametros(2));

% resultados para los 3 ejemplos:
% ejemplo 2
% ----
% funcion objetivo = -1.899024
% x1= 0.975841 - x2= 1.006707
% deltf (10 mejores): 0.0424034
% ----
% funcion objetivo = -1.981827
% x1= 0.990708 - x2= 0.994962
% deltf (10 mejores): 2.22045e-15
% ----
% funcion objetivo = -1.995138
% x1= 0.995446 - x2= 1.003017
% deltf (10 mejores): 0.140318
%------------------------------
% ejemplo 3
% ----
% funcion objetivo = 2.000028
% x1= -1.563292 - x2= 0.436736
% deltf (10 mejores): 3.85901e-05
% ----
% funcion objetivo = 2.000015
% x1= -1.810381 - x2= 0.189605
% deltf (10 mejores): 4.02508e-05
% ----
% funcion objetivo = 2.000000
% x1= -1.307977 - x2= 0.692023
% deltf (10 mejores): 5.84853e-05
%------------------------------
% ejemplo 1
% ----
% funcion objetivo = 0.040000
% x1= 0.800000 - x2= 0.639928
% deltf (10 mejores): 1.2454e-05
% ----
% funcion objetivo = 0.040000
% x1= 0.800000 - x2= 0.639879
% deltf (10 mejores): 9.0546e-06
% ----
% funcion objetivo = 0.040006
% x1= 0.799986 - x2= 0.639731
% deltf (10 mejores): 1.74389e-05
