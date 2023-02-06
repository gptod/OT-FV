clear
close all


% set test case fixing intitial and final density.
% See function boundary_bc for the available options.
test_cases = ["sin"];%,"compression"];
for  test_case = test_cases;
	test_case = test_case;
	disp(test_case)

	% select the directory where to store the results
	folder_runs='runs_BenamouBrenier';
	if ~isfolder(folder_runs)
		mkdir(folder_runs);
	end

	%
	% set some globals controls for reproducing experiments
	%
	plot_figures = 1; % print rho evolutions
	compute_err = 1; % compute errors with respect to exact solutions


	%
	% SPATIAL DISCRETIZATION
	%

	
	% refine level. Available from 1 to 5
	for h_i = 1;

		% recostruction used
		% rec == 1 : linear
		% rec == 2 : harmonic (some prec.s do not support it)
		rec = 1;

		% use a refined grid for phi
		% 0 : one grid are used
		% 1 : tow grids are used
		twogrids = 1;

		% create folder storing data for test case
		folder=sprintf('%s/%s_rec%d/',folder_runs,test_case,rec);
		if ~isfolder(folder)
			mkdir(folder);
		end
		
		
		% grids for rho and phi (see TPFA_grid class)
		% are read from file
		
		fn1 = 'meshes/mesh_gen/tri2/mesh1/coord.txt'; %coordinates
		fn2 = 'meshes/mesh_gen/tri2/mesh1/topol.txt'; %topology
		[nodes,cells] = get_mesh(fn1,fn2,0,0,0,0);
		grid_rho = TPFA_grid;
		grid_rho.init(cells,nodes);
		
		%
		% TEMPORAL DISCRETIZATION delta=1/N
		%
		for dt_i = 1
			%for dt_i = h_i +1
			% number of time steps
			%Ntime = 4*(2^(dt_i));
            Ntime = 4*(2^(dt_i))-1;

			
			%
			% PROBLEM DATA
			%
			% Set initial and final rho density
			% In same cases the exact solution is known
			[rho_in,rho_f,...
			 mass,midpoint,exact_phi_function,exact_rho_function,exact_W2,bc_sol] = ...
			bc_density(test_case,grid_rho.cc,grid_rho.area);

			
			
			% solve with interior point
			[phi,rho,ierr,info_solver] = solve_benamoubrenier(cells,nodes,Ntime,...
															  rho_in,rho_f,...
															  1, rec, 1e-5);
			
			% plot
			if (plot_figures)
                Ndsp=3;
                dsp=ceil(Ntime/(Ndsp+1)):ceil(Ntime/(Ndsp+1)):ceil(Ntime-Ntime/(Ndsp+1));
				plot_rhos_2d(grid_rho,rho_in,rho,rho_f,dsp)
			end
		end % mesh size
	end % time size
end  % test case


