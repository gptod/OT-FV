function err_rhoth = compute_error_rho(rho,rho_in, rho_final, grid_rho,exact_rho)
	%
	Nr= size(rho_in,1);
	N = size(rho,1)/Nr;
	
	% Compute the errors
  err_rhoth = 0;
  for k=1:N
    % error on the geodesic
    err_rhoth = err_rhoth + (1/(N+1))*sum(grid_rho.area.*abs(rho(1+(k-1)*Nr:k*Nr)-exact_rho(1+(k-1)*Nr:k*Nr)));
  end
end
