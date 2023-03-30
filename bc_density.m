function [rho_in,rho_f,mass,midpoint,potential,geodesic,W2,bc_sol] = bc_density(test_case,cc,area)

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
  
  midpoint = [];
  potential = [];
  geodesic = [];
  W2 = [];
  bc_sol = 0;

elseif (strcmp(test_case,'gauss_plus'))
    
	% Gaussian densities
	lift = 1e-3;
  initial =@(x,y) exp((-(x-0.3).^2-(y-0.3).^2)/0.05)+lift;
  final =@(x,y) exp((-(x-0.7).^2-(y-0.7).^2)/0.05)+lift;
  rho_in = initial(cc(:,1),cc(:,2));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2));
  rho_f = rho_f*mass/sum(rho_f.*area);
  
  midpoint = [];
  potential = [];
  geodesic = [];
  W2 = [];
  bc_sol = 0;
elseif (strcmp(test_case,'gauss_wide'))
    
	% Gaussian densities
	std_variation = 0.1;
  initial =@(x,y) exp((-(x-0.3).^2-(y-0.3).^2)/std_variation);
  final =@(x,y) exp((-(x-0.7).^2-(y-0.7).^2)/std_variation);
  rho_in = initial(cc(:,1),cc(:,2));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2));
  rho_f = rho_f*mass/sum(rho_f.*area);
  midpoint = [];
  potential = [];
  geodesic = [];
  W2 = [];
  bc_sol = 0;
	
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
    
  midpoint = @(x,y) ((sqrt((x-0.5).^2+(y-0.5).^2)-0.3)<=0).*(1+cos(supp*((x-0.5).^2+(y-0.5).^2)));
  potential=@(x,y,t) 0.4*x+0.4*y-(0.4^2)*t-0.4;
  geodesic =@(x,y,t) ((sqrt((x-0.3-t*0.4).^2+(y-0.3-t*0.4).^2)-0.3)<=0).*(1+cos(supp*((x-0.3-t*0.4).^2+(y-0.3-t*0.4).^2)));
  W2 = sqrt(2)*0.4;
  bc_sol = 1;
  
elseif (strcmp(test_case,'sin_old'))
    
% Sinusoidal functions
  initial =@(x,y) cos(2*pi*sqrt((x-0.5).^2+(y-0.5).^2))+1.5; 
  final =@(x,y) -cos(2*pi*sqrt((x-0.5).^2+(y-0.5).^2))+1.5;
  rho_in = initial(cc(:,1),cc(:,2));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2));
  rho_f = rho_f*mass/sum(rho_f.*area);
  
  midpoint = [];
  potential = [];
  geodesic = [];
  W2 = [];
  bc_sol = 0;

elseif (strcmp(test_case,'cross'))
				% % Cross distributed densities
  beta = pi/4; R = [cos(beta) -sin(beta); sin(beta) cos(beta)];
  ccrot = ([cc(:,1) cc(:,2)]-0.5)*R+0.5;
  cross =@(x,y) (x>=0.4).*(x<=0.6).*(y>=0).*(y<=1)+(y>=0.4).*(y<=0.6).*(x>=0).*(x<=1);
  rho_in = cross(cc(:,1),cc(:,2));
  rho_in(rho_in>0) = rho_in(rho_in>0)./rho_in(rho_in>0);
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = cross(ccrot(:,1),ccrot(:,2));
  rho_f(rho_f>0) = rho_f(rho_f>0)./rho_f(rho_f>0);
  rho_f = rho_f*mass/sum(rho_f.*area);

  midpoint = [];
  potential = [];
  geodesic = [];
  W2 = [];
  bc_sol = 0;

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
    
  midpoint =@(x,y) ((((x-0.5)/(0.5*sc+0.5)).^2-0.5^2)<=0).*(1+cos(2*pi*((x-0.5)/(0.5*sc+0.5))))/(0.5*sc+0.5);
  potential =@(x,y,t)((sc-1)/((sc-1)*t+1))*0.5*(x-0.5).^2;
  geodesic =@(x,y,t) ((abs(x-0.5)-0.5*(t*(sc-1)+1))<=0).*(1+cos(2*pi*((x-0.5)/(t*(sc-1)+1))))/(t*(sc-1)+1);
  W2 = sqrt(((pi^2-6)*(sc-1)^2)/(12*pi^2));
  bc_sol = 1;
  
elseif (strcmp(test_case,'gauss_3d'))
  
  % Gaussian densities
  initial =@(x,y,z) exp((-(x-0.3).^2-(y-0.3).^2-(z-0.3).^2)/0.1);
  final =@(x,y,z) exp((-(x-0.7).^2-(y-0.7).^2-(z-0.7).^2)/0.1);
  rho_in = initial(cc(:,1),cc(:,2),cc(:,3));
  mass = sum(area.*rho_in);
  rho_in = rho_in/mass; mass = sum(area.*rho_in); %normalization
  rho_f = final(cc(:,1),cc(:,2),cc(:,3));
  rho_f = rho_f*mass/sum(rho_f.*area);
  
  midpoint = [];
  potential = [];
  geodesic = [];
  W2 = [];
  bc_sol = 0;
  
else
  disp('Testcase not supported')
end
