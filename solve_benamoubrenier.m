function [phi,rho,ierr,info_solver] = solve_benamoubrenier(topol,coord,Ntime,...
														   rho_initial,rho_final,...
														   twogrids, rec, tol,...
														   phi_initial_guess,rho_initial_guess)

	% assemble structures with geometrical info 
	grid_rho = TPFA_grid;
	grid_rho.init(topol,coord);
	if twogrids
		[grid_phi, I] = grid_rho.refine();
	else
		grid_phi = grid_rho;
		I = speye(grid_rho.ncell);
	end

	% set controls for Interior Point solver and linear solver
	IP_ctrl=IP_controls;
	IP_ctrl.eps_0 = tol;
	IP_verbose = 2;
	[linear_solver_ctrl,label] = set_linear_algebra_ctrl('direct',rec);

	
	% check data and set initial solution
	mass = grid_rho.area'*rho_initial;
	mass_final = grid_rho.area'*rho_final;
	
	if abs( mass - mass_final) > 1e-12
		error('Initial and final denisity have different mass')
	end
	
	Nr = grid_rho.ncell * Ntime;
	Np = grid_phi.ncell * (Ntime + 1);
	
	if exist('phi_initial_guess ','var') 
		phi0 = phi_initial_guess;
	else
		phi0 = ones(Np,1);
	end

	if exist('rho_initial_guess ','var') 
		rho0 = rho_initial_guess;
	else
		rho0  = ones(Nr,1);
	end
	s0 = 1.0./rho0;
	
	phi_rho_slack_initial_guess = [phi0;rho0;s0];
		

	% run solver
	
	[phi,rho,slack,W2th,info_solver] = l2otp_solve(grid_rho, grid_phi,I, rec, Ntime,...
												   IP_ctrl,linear_solver_ctrl,... 
												   rho_initial,rho_final, phi_rho_slack_initial_guess);

	ierr = info_solver.ierr;

end
