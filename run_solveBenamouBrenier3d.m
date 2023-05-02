clear
close all

test_case = "gauss_3d";

% select the directory where to store the results
folder_runs='runs_BenamouBrenier3d';
if ~isfolder(folder_runs)
	mkdir(folder_runs);
end


%
% set some globals controls for reproducing experiments
%
plot_figures = 1; % print rho evolutions
compute_err = 0; % compute errors with respect to exact solutions


%
% SPATIAL DISCRETIZATION
%


% mesh level. Available from 0 to 2
for h_i = 0
h = 2^(3+h_i);
% TEMPORAL DISCRETIZATION delta=1/(N+1)
dt = h_i;
Ntime = 2^(3+h_i)-1;

% recostruction used
% rec == 1 : linear
% rec == 2 : harmonic (some prec.s do not support it)
rec = 1;

% use a refined grid for phi
% 0 : one grid are used
% 1 : tow grids are used
twogrids = 0;

% create folder storing data for test case
folder=sprintf('%s/%s_rec%d/',folder_runs,test_case,rec);
if ~isfolder(folder)
	mkdir(folder);
end

test_case_label = sprintf('%s_rec%d_h%d_N%d_',test_case,rec,h_i,dt);

% grids for rho and phi (see TPFA_grid class)
% are read from file
[grid_rho, grid_phi, I] = init_grids_3d(1, h);
		
			
%
% PROBLEM DATA
%
% Set initial and final rho density
% In same cases the exact solution is known
[rho_in,rho_f,...
 mass,midpoint,exact_phi_function,exact_rho_function,exact_W2,bc_sol] = ...
bc_density(test_case,grid_rho.cc,grid_rho.area);

			

% set controls for Interior Point solver and linear solver
IP_ctrl=IP_controls;



IP_ctrl.eps_0 = 1e-5;
IP_verbose = 2;
approach='primal';
[linear_solver_ctrl,labels] = set_linear_algebra_ctrl(approach,rec);
folder_approach=strcat(folder,'/',approach,'/');
if ~isfolder(folder_approach)
	mkdir(folder_approach);
end


linear_algebra_label=labels(1,1);
experiment_label = strcat(test_case_label,linear_algebra_label)
filename = strcat(folder_approach,experiment_label)
% csv file with performance resume

IP_ctrl.save_csv=1;
IP_ctrl.file_csv=strcat(filename,'.csv');

IP_ctrl.save_log=1;
IP_ctrl.file_log=strcat(filename,'.log');
	
Nr = grid_rho.ncells * Ntime
Np = grid_phi.ncells * (Ntime + 1)
	
phi0 = ones(Np,1);
rho0  = ones(Nr,1);
s0 = 1.0./rho0;
phi_rho_slack_initial_guess = [phi0;rho0;s0];
		

% run solver
[phi,rho,slack,W2th,info_solver] = l2otp_solve(grid_rho, grid_phi,I, rec, Ntime,...
												   IP_ctrl,linear_solver_ctrl,... 
												   rho_in,rho_f, phi_rho_slack_initial_guess);

ierr = info_solver.ierr;



			
% plot
if (plot_figures)
    Ndsp=3;
    dsp=ceil(Ntime/(Ndsp+1)):ceil(Ntime/(Ndsp+1)):ceil(Ntime-Ntime/(Ndsp+1));
	plot_rhos(grid_rho,rho_in,rho,rho_f,dsp)
end


end 
