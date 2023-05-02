function [err_cost,err_p,err_rhoth]=compute_errors(W2,rho_in,rho,rho_final,phi,... % variable 
																									 grid_rho,grid_phi,RHt,It,... % geometry 
																									 mass,exact_W2,exact_rho,exact_phi) % exact formulas
	%
	% compute errors w.r.t to exact solution given
	% exact W2 distance, formulas for phi and rho.
	%  
	
	indc=1; % fix this
  % shift the constant of phi
  [~,indc]=min((grid_phi.cc(:,1)-0.5).^2+(grid_phi.cc(:,2)-0.5).^2);
  phi = phi+(exact_phi(0.5,0.5,0)-phi(indc));

	rho_all=[rho_in;rho;rho_final];
	rhoa = RHt*It*rho_all;
	
	N=size(rho,1)/grid_rho.ncells;
	Nt=N+1;
	
  % Compute the errors
  err_cost = abs(W2-exact_W2); % cost error
  err_p = 0;
  err_rhoth = 0;
  for k=1:Nt
    t = (k-0.5)/Nt;
    phi_real = exact_phi(grid_phi.cc(:,1),grid_phi.cc(:,2),t);
    rho_real = exact_rho(grid_rho.cc(:,1),grid_rho.cc(:,2),t);
    rho_real = rho_real*mass/sum(rho_real.*grid_rho.area);
    rho_t = 0.5*(rho_all((k-1)*grid_rho.ncells+1:k*grid_rho.ncells)+rho_all(k*grid_rho.ncells+1:(k+1)*grid_rho.ncells));
    % error on the potential
    err_p = err_p + (1/Nt)*sum( rhoa((k-1)*grid_phi.ncells+1:k*grid_phi.ncells).*grid_phi.area.*(phi_real-phi((k-1)*grid_phi.ncells+1:k*grid_phi.ncells)).^2 );
    % error on the geodesic
    err_rhoth = err_rhoth + (1/Nt)*sum(grid_rho.area.*abs(rho_t-rho_real));
  end
  err_p = sqrt(err_p);
end
