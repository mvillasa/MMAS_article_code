% author: Saul Buitrago (USB)
function eval_funcion = fun_objetive_2p(z,example)
% evalua la funcion objetivo (error) en z (2 parametros)
switch example
  case 3
    %---- funcion3 en [-2,0]x[0,2]
    eval_funcion = abs(-z(1)+z(2)-2.)+2.;
  case 2
    %---- funcion2 en [0,2]x[0,2]
    eval_funcion = (z(1)-1.)^2+(z(2)-1.)^2-cos(18*(z(2)-1.))-cos(18*(z(1)-1.));
  case 1
    %----- funcion1 [-1.2,0.8]x[-1.2,0.8]
    eval_funcion = 10.*(z(2)-z(1)*z(1))^2 + (1.-z(1))^2;
end
end

