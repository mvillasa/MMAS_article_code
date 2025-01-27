% author: Saul Buitrago (USB)
function [xx,fxy,deltf10,varianza,icont5] = ...
       opt_ex_in(poblacion,input_no,cycle,example)
% minimizador del error
fprintf('poblacion: %d - param: %d - ciclos: %d\n',poblacion,input_no,cycle);
 
tol = 1.e-4; %1e-8; % tolerance for deltf  1.e-3; 1.e-6
alpha = 1.2; % parametro para estirar el radio "irad" de las bolas
fact_red = 0.9; %sebb 0.5; %0.6; %0.2 0.4  0.6  0.8  1.0;
nit_rate = 0.25; % rate to generate new individuals per step 0.25
nit = floor(poblacion*nit_rate); % numero de individuos por iteracion
fprintf('numero de individuos por ciclo: %d\n',nit);
pendie = 0.005; % tolerance slope of the plane
% Define function f and fun_obj evaluation of individuals
switch input_no
  case 2
    nmaxp = 500000; %500000; % base para seleccionar en maximo esparcido (maximum spread)
    switch example
        case 3
         lb=[-2.0,0.0]; ub=[0.0,2.0]; % funcion 3  [-2,0]x[0,2]
        case 2
         lb=[0.0,0.0]; ub=[2.0,2.0]; % funcion 2  [0,2]x[0,2]
        case 1
         lb=[-1.2,-1.2]; ub=[0.8,0.8]; % funcion 1  [-1.2,0.8]x[-1.2,0.8]
    end
    f = @fun_objetive_2p;
  otherwise
    disp('number of variables only 2');
    return
end
% inicializacion de arreglos
xm = zeros(poblacion+(cycle+1)*nit,input_no);
f_o = zeros(1,poblacion+(cycle+1)*nit);
error1 = zeros(1,cycle); %error2 = zeros(1,cycle);
fprintf('\n\n**initial population**\n');

% Populate the initial population by maximum spread 
%   (generar la población inicial por dispersión máxima)
ix = zeros(nmaxp,input_no);
for i = 1:nmaxp
    for j =1:input_no
            aux1 = ub(j);
            if aux1 >= lb(j)
                ix(i,j) = (aux1 - lb(j))*rand(1,1) + lb(j);
            else  % caso cuando no hay suficiente gas disponible
                ix(i,j) = 0.;
            end
    end
end
% calculo del radio xrad de la esfera
xnt = 0.;
for i=1:input_no
    xnt = xnt + log(ub(i)-lb(i));
end
xnv = 1./input_no;
xnt = (xnt - log(poblacion))*xnv;
xrad = exp(xnt)/2.;
xrad_dil = xrad*alpha;
% generacion de los "poblacion" puntos a partir de los nmaxp puntos
%   verificando que esten fuera de las esferas de radio xrad_dil y
%   centro los puntos generados hasta el momento
item_aux = randi([1,nmaxp],1,nmaxp);
ite = item_aux(1);
xm(1,:) = ix(ite,:);
f_o(1) = f(xm(1,:),example); % llamado a f
param_fact_red = 1;
while param_fact_red == 1
    param_fact_red = 0;
    xrad_dil = xrad_dil*fact_red;
    icont = 1;
    for i=2:poblacion
        icont = icont +1;
        if icont >= nmaxp
            fprintf('Se alcanzo el maximo de %d, Aumentar el valor de nmaxp, caso 1\n', nmaxp);
            %param_fact_red = 1;
            return;
        end
        inorm = 0.;
        while inorm < xrad_dil
           ite = item_aux(icont);
           item = ix(ite,:);
           for j=1:i-1
               inorm = norm(item-xm(j,:),inf);
               if inorm < xrad_dil
                   break;
               end
           end
           icont = icont + 1;
           if icont >= nmaxp
               fprintf('Se alcanzo el maximo de %d, Aumentar el valor de nmaxp, caso 2\n', nmaxp);
               return;
           end
        end
        if icont >= nmaxp
            fprintf('Se alcanzo el maximo de %d, Aumentar el valor de nmaxp, caso 2\n', nmaxp);
            return;
        else
            xm(i,:) = item;
            f_o(i) = f(xm(i,:),example); % llamado a f
        end 
    end
end  %linea 1054
fprintf('maximum spread\n');
for j=1:input_no
    fprintf('coord %d: ',j);
    for ii=1:poblacion, fprintf('%f ',xm(ii,j)); end
    fprintf('\n');
end
fprintf('f_o:    ');
for ii=1:poblacion, fprintf('%f ',f_o(ii)); end
fprintf('\n');
if input_no == 2    % radio esfera xrad_dil
    figure(2);
    plot(xm(1:poblacion,1),xm(1:poblacion,2),'r*');
    axis([lb(1), ub(1), lb(2), ub(2)]);
    hold on;
    % dibujar las esferas de centro xm y radio xrad_dil
    centers = xm(1:poblacion,:);
    radii = ones(poblacion,1)*xrad_dil;
    viscircles(centers,radii,'EdgeColor','b','LineStyle',':','LineWidth',0.1);
    %axis([0 gas_constraint 0 gas_constraint])
    %voronoi(xm(1:poblacion,1),xm(1:poblacion,2));
    %fill([1000,1000,0],[0,1000,1000],[1,1,1]);
    title('maximun spread');
    hold off;
end
% se ordena el vector inicial f_o de menor a mayor
[~,index] = sort(f_o(1:poblacion),'ascend');
imax = index(1); fmax = f_o(imax);
imin = index(poblacion); fmin = f_o(imin); 
fprintf('min y max: %d %f ** %d %f\n',imax,fmax,imin,fmin);
%fprintf('press Intro to continue\n'); pause;

% calculo de la mayor diferencia de f_o
deltf = abs(fmax - fmin); disp(deltf);
ip = 1;
%icont2 = 0;
icont5 = poblacion; % fin calculo de puntos iniciales
err1 = 10000;
fun_obj = zeros(1,cycle); no = zeros(1,cycle);
%pause;
while (ip <= cycle) && (deltf > tol)
    [~,indd] = sort(f_o(1:icont5),'ascend');
    aux_f_o = f_o(indd(1:poblacion));
    aux_xm = xm(indd(1:poblacion),:);
    no(ip) = ip; fun_obj(ip) = f_o(indd(1)); % se almacenan el error por iteracion
    % se imprimen los 20 mejores
    disp('los 20 mejores: aux_xm'); fprintf('inicio ip=%d\n',ip);
    for j=1:input_no
      fprintf('coord %d: ',j);
      for ii=1:poblacion, fprintf('%f ',aux_xm(ii,j)); end
      fprintf('\n');
    end
    fprintf('f_o:    ');
    for ii=1:poblacion, fprintf('%f ',aux_f_o(ii)); end
    fprintf('\n'); %pause;
    %disp('aux_xm'); disp([aux_xm,aux_f_o']); pause;
    rs = factorial(poblacion)/(factorial(poblacion-3)*3*2); % numero de combinaciones de n elementos tomadas de 3 en 3 (n!/((n-3)! * 3!))
    rs = ceil(rs);
    item_aux1 = randperm(rs); % row vector containing a random permutation of the integers from 1 to rs inclusive
    contador = 1;
    array_comb = zeros(1000,3); % linea 1321
    for i1=1:poblacion
       for i2=i1+1:poblacion
          for i3=i2+1:poblacion
              array_comb(contador,1) = i1;
              array_comb(contador,2) = i2;
              array_comb(contador,3) = i3;
              % fprintf('%d %d %d \n',i1,i2,i3);
              if contador > rs  % nunca debe suceder
                 fprintf('---- contador %d > rs %d\n',rs, contador);
                 return;
              end
              contador = contador + 1;
          end
       end
    end
    % rs debe ser contador-1 siempre
    % fprintf('rs %d ** contador-1 %d\n',rs,contador-1);
%     for i1=1:contador-1
%        fprintf('%d %d %d\n',array_comb(i1,1:3)); 
%     end
    icont1 = 0; % icont1: number of points generated per step
    icont2 = 1; % icont2 - 1: number of function evaluation carried out per step

    while icont2 <= nit % numero de individuos por iteracion
      genera = 1;
      while genera
        icont1 = icont1 + 1;
        if icont1 > rs
          fprintf('se agoto el conjunto de combinaciones para generar nuevos puntos\n');
          fprintf('ip: %d - icont1: %d - combinaciones %d \n',ip,icont1,contador-1);
          %icont1 = 0;
          fprintf('press Intro to continue\n'); %pause; 
          return; %continue; %SEBB202312
        end
        % seleccion de 3 puntos del conjunto de "poblacion" de puntos
        disp('seleccion de 3 puntos');
        item = item_aux1(icont1);
        itemp = zeros(3,input_no);
        ftemp = zeros(1,3);
        for i=1:3    % linea 1345 
            aux_int = array_comb(item,i);
            itemp(i,1:input_no) = aux_xm(aux_int,1:input_no); % aux_xm tiene los "poblacion" mejores
            fprintf('%d ** ',aux_int);
            for ii=1:input_no
                fprintf('%f ',itemp(i,ii));
            end
            ftemp(i) = aux_f_o(aux_int);
            fprintf('** %.12f\n',ftemp(i));
        end % linea 1355
        [~,iauxx] = sort(ftemp,'ascend');
        aux1 = norm(itemp(iauxx(1),1:input_no)-itemp(iauxx(3),1:input_no),2); %disp(aux1);
        aux2 = abs(ftemp(iauxx(1)) - ftemp(iauxx(3))); %disp(aux2);
        pend = aux2/aux1;
        if pend == 0
           fprintf('** ip: %d - pend es cero',ip); pause;
        end
        fprintf('pend: %f\n',pend);
        inew = sum(itemp(iauxx(1:2),1:input_no))/2; disp(inew); % (x1+x2)/2
        dentro = 1;
        % generacion del nuevo punto inew
        centro = 0; aux = 0;
        if pend < pendie  % linea 1376
            % centro de xg y de x3 con xg=(x1+x2)/2
            inew = (inew + itemp(iauxx(3),1:input_no))/2;
            centro = 1; ii=0;
        else
            % imagen de x3 con respecto a xg, verificando que esta en el dominio % linea 1383
            dentro = 0; nn = 25; refl = linspace(2,1.1,nn); %SEBB nn=25
            for ii=1:nn
               inew1 = itemp(iauxx(3),1:input_no) + refl(ii)*(inew - itemp(iauxx(3),1:input_no));
               aux = all(and(inew1>lb,inew1<ub)); fprintf('aux: %d\n',aux);
               if aux == 1
                   dentro = 1; centro = 0;
                   break;
               end
            end
            inew = inew1;
            fprintf('ii: %d - aux: %d\n',ii,aux); disp(inew);
            if ii == nn && aux == 0 % ii alcanzo el maximo e inew1 no esta en [0,ub]
                disp('se alcanzo el maximo del lazo en: imagen de x3 con respecto a xg');
                inew1 = itemp(iauxx(3),1:input_no) + 0.85*(inew - itemp(iauxx(3),1:input_no));
                aux = all(and(inew1>=0,inew1<ub)); %si aux=1 entonces inew1 esta en [0,ub]
                %fprintf('aux: %d\n',aux);
                inew = inew1;
                % se genero un inew dentro del triangulo (factor 0.85, 0.95)
                %pause;
            end
        end
        disp('-- nuevo generado'); 
        fprintf('ii: %d - dentro: %d - centro %d - aux: %d\n',ii,dentro,centro,aux); 
        disp(inew);
        if dentro == 0
            genera = 1;
            continue
        end
        disp('-- nuevo generado'); disp(inew);
        genera = 0;
        if input_no == 2
          figure(3); % se grafica la generacion de los puntos nuevos
          plot(aux_xm(1:poblacion,1),aux_xm(1:poblacion,2),'r*');
          axis([lb(1),ub(1),lb(2),ub(2)]);  %sebb era [0,ub(1),0,ub(2)]
          hold on;
          plot(itemp(iauxx(1:2),1),itemp(iauxx(1:2),2),'-k');
          plot(itemp(iauxx(3),1),itemp(iauxx(3),2),'ok');
          plot(inew(1),inew(2),'ob');
          hold off;
        end
        xm(poblacion+(ip-1)*nit+icont2,:) = inew;
        % evaluacion en f de los nuevos puntos inew
        f_o(poblacion+(ip-1)*nit+icont2) = f(inew,example); % llamado a f
      end % fin while genera

      fprintf('**** icont2: %d - nit: %d - ip: %d\n',icont2,nit,ip);
      icont2 = icont2 + 1; 
      % se debe ordenar xm en funcion de f_o para trabajar con los
      % "poblacion" mejores en la iteracion siguiente
      % crear un arreglo paralelo aux_xm a xm con los "poblacion" mejores
      % se hace al principio del "while (ip <= cycle) && (deltf > tol)"
      
    end % fin del while icont2 <= nit
    icont5 = icont5 + nit;
    %[~,indd] = sort(f_o(1:poblacion+(ip-1)*nit+icont2),'ascend');
    [~,indd] = sort(f_o(1:icont5),'ascend');
    %maxf_o = f_o(indd(1)); minf_o = f_o(indd(poblacion));
    minf_o = f_o(indd(1)); maxf_o = f_o(indd(poblacion)); %para minimizar f
    err1 = minf_o;
    error1(ip) = abs(maxf_o-minf_o); deltf = error1(ip);
    fprintf('**** ip %d - deltf %g\n',ip,deltf);
    fprintf('error de los %d mejores: \n',poblacion);
    if input_no ~= 2
      fprintf('%g ',f_o(indd(1:poblacion)));
      fprintf('\n');
    else
      for iii=1:poblacion
        fprintf('%f %f %g\n',xm(indd(iii),1),xm(indd(iii),2),f_o(indd(iii))); 
      end
    end
    fprintf('\npress Intro to continue\n'); %pause;
    ip = ip + 1;
end % fin del while (ip<=cycle) && (err1>tol) && (deltf>tol)
n_ite = ip-1;

fprintf('ip: %d - err1: %g - deltf: %g\n',ip,err1,deltf);
fprintf('\n****RESULTADOS****\n');
[minimof,iminf] = min(f_o(1:icont5));
fprintf('%d parametros: ',iminf);
for i=1:input_no
    fprintf('%.15f ',xm(iminf,i));
end
xx = xm(iminf,1:input_no);
fxy = minimof;
fprintf('\n');
fprintf('%d error: %g\n',iminf,minimof);
fprintf('deltf (poblacion): %g\n',deltf);
[~,index] = sort(f_o(1:icont5),'ascend'); 
fmin = f_o(index(1));
fmax = f_o(index(10)); 
deltf10 = abs(fmax - fmin); % mayor diferencia de los 10 mejores f_o
fprintf('deltf (10 mejores): %g\n',deltf10);
varianza = var(f_o(index(1:10)));
fprintf('varianza: %g\n',varianza); % en notacion exponencial usar %e
fprintf('function evaluations: %d\n',icont5);
fprintf('numero de generaciones: %d\n',n_ite);
if input_no == 2  %sebb solo para caso de 2 parametros
    figure(4); % se grafican todos los puntos
    plot(xm(1:icont5,1),xm(1:icont5,2),'r*');
    axis([lb(1),ub(1),lb(2),ub(2)]);
    hold on;
    plot(xm(index(1:10),1),xm(index(1:10),2),'bo');
    % ploting the optimum
    x = xm(iminf,1); y = xm(iminf,2);
    plot(x,y,'k*');
    hold off;
end
figure(5);
plot(no(1:n_ite),fun_obj(1:n_ite),'Linewidth',2,'Color','r');
title('evolucion de la funcion objetivo','FontSize',10);
ylabel('fun-obj','FontSize',10);
xlabel('Iteration','FontSize',10);

end