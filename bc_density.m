function [rho_in,rho_f,mass] = bc_density(test_case,cc,area)

% Function that generates the discrete boundary conditions on the density

if (strcmp(test_case,'gauss'))
				% Gaussian densities
  initial =@(x,y) exp((-(x-0.3).^2-(y-0.3).^2)/0.05);
  final =@(x,y) exp((-(x-0.7).^2-(y-0.7).^2)/0.05);
  rho_in = initial(cc(:,1),cc(:,2));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2));
  rho_f = rho_f*mass/sum(rho_f.*area);
elseif (strcmp(test_case,'sin'))
  % convergence test 1
  supp = pi/(0.3^2);
  initial =@(x,y) ((sqrt((x-0.3).^2+(y-0.3).^2)-0.3)<=0).*(1+cos(supp*((x-0.3).^2+(y-0.3).^2)));
  final =@(x,y) ((sqrt((x-0.7).^2+(y-0.7).^2)-0.3)<=0).*(1+cos(supp*((x-0.7).^2+(y-0.7).^2)));
  %supp = pi/0.3;
  %initial =@(x,y) (((x-0.3).^2+(y-0.3).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.3).^2+(y-0.3).^2))));
  %final =@(x,y) (((x-0.7).^2+(y-0.7).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.7).^2+(y-0.7).^2))));
  rho_in = initial(cc(:,1),cc(:,2));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2));
  rho_f = rho_f*mass/sum(rho_f.*area);
elseif (strcmp(test_case,'sin_old'))
% Sinusoidal functions
  initial =@(x,y) cos(2*pi*sqrt((x-0.5).^2+(y-0.5).^2))+1.5; 
  final =@(x,y) -cos(2*pi*sqrt((x-0.5).^2+(y-0.5).^2))+1.5;
  rho_in = initial(cc(:,1),cc(:,2));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2));
  rho_f = rho_f*mass/sum(rho_f.*area);

elseif (strcmp(test_case,'cross'))
				% % Cross distributed densities
  beta = pi/4; R = [cos(beta) -sin(beta); sin(beta) cos(beta)];
  ccrot = ([cc2h(:,1) cc2h(:,2)]-0.5)*R+0.5;
  cross =@(x,y) (x>=0.4).*(x<=0.6).*(y>=0).*(y<=1)+(y>=0.4).*(y<=0.6).*(x>=0).*(x<=1);
  rho_in = cross(cc2h(:,1),cc2h(:,2));
  rho_in(rho_in>0) = rho_in(rho_in>0)./rho_in(rho_in>0);
  mass = sum(area2h.*rho_in);
  rho_in = rho_in/mass; mass = sum(area2h.*rho_in); %normalization
  rho_f = cross(ccrot(:,1),ccrot(:,2));
  rho_f(rho_f>0) = rho_f(rho_f>0)./rho_f(rho_f>0);
  rho_f = rho_f*mass/sum(rho_f.*area2h);


elseif( strcmp(test_case,'compression'))
    % convergence test 2
  sc = 0.3;
  initial =@(x,y) (1+cos(2*pi*(x-0.5)));
  final =@(x,y) ((((x-0.5)/sc).^2-0.5^2)<=0).*(1+cos(2*pi*((x-0.5)/sc)))/sc;
  rho_in = initial(cc(:,1),cc(:,2));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2));
  rho_f = rho_f*mass/sum(rho_f.*area);
  
else
  disp('Testcase not supported')
end
