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
% start the solution from prebuilt initial guess
% 
read_from_file = 0; % 0 : standard initial guess >0: start from the corresponding IP iterations
folder_restart='initial_solution'; % location of initial solution. See more below

%
% set some globals controls for reproducing experiments
%
plot_figures = 1; % print rho evolutions
compute_err = 1; % compute errors with respect to exact solutions


%
% SPATIAL DISCRETIZATION
%

% Different type of meshes
% 1: 
% 2:
% 3:
% 4,5,6: like 1,2,3 but with two grid level
for mesh_type = 5;
	
	% refine level. Available from 1 to 5
	for h_i = 3;

		% recostruction used
		% rec == 1 : linear
		% rec == 2 : harmonic (some prec.s do not support it)
		rec = 1;

		% create folder storing data for test case
		folder=sprintf('%s/%s_rec%d/',folder_runs,test_case,rec);
		if ~isfolder(folder)
			mkdir(folder);
		end
		folder=sprintf('%s/%s_rec%d/initial_solution/',folder_runs,test_case,rec);
		if ~isfolder(folder)
			mkdir(folder);
		end
		
		
		% grids for rho and phi (see TPFA_grid class)
		% are read from file
		twogrids = 1;
		fn1 = 'meshes/mesh_gen/tri2/mesh1/coord.txt'; %coordinates
		fn2 = 'meshes/mesh_gen/tri2/mesh1/topol.txt'; %topology
		[nodes,cells] = get_mesh(fn1,fn2,0,0,0,0);
		grid_rho = TPFA_grid;
		grid_rho.init(cells,nodes);
		if twogrids
			disp('here')
			%grid_phi = TPFA_grid;
			[grid_phi, I] = grid_rho.refine();
			disp('refined')
		else
			grid_phi = grid_rho;
			I = speye(grid_rho.ncell);
		end
	
		%
		% TEMPORAL DISCRETIZATION delta=1/N
		%
		for dt_i = 1
		%for dt_i = h_i +1
			% number of time steps
			Ntime = 4*(2^(dt_i));

			% set problem dimension
			ncell_phi = grid_phi.ncell;
			ncell_rho = grid_rho.ncell;
			Np = (Ntime + 1) * ncell_phi;
			Nr = (Ntime) * ncell_rho;
			


			%
			% PROBLEM DATA
			%
			% Set initial and final rho density
			% In same cases the exact solution is known
			[rho_in,rho_f,...
			 mass,midpoint,exact_phi_function,exact_rho_function,exact_W2,bc_sol] = ...
			bc_density(test_case,grid_rho.cc,grid_rho.area);

			% store exact solution 
			phi_rho_slack_reference_solution = zeros(Np+2*Nr,1);
			%exact_rho_vec=zeros(Nr,1);
			if ( bc_sol == 1)
				for k=1:Ntime
					t=k/(Ntime+1);
					rho_real = exact_rho_function(grid_rho.cc(:,1),grid_rho.cc(:,2),t);
					rho_real = rho_real*mass/sum(rho_real.*grid_rho.area);
					phi_rho_slack_reference_solution(Np+1+(k-1)*ncell_rho:Np+k*ncell_rho)=rho_real;
				end
			end			
		

			% set test case label
			test_case_label = sprintf('%s_rec%d_mesh%d_h%d_N%0.5d_',test_case,rec,mesh_type,h_i,Ntime);


			%
			% ALGORITHM CONTROLS
			%
			% Init Interior Point controls.
			% Defalt values are set in IP_controls class
			IP_ctrl=IP_controls;
			
			% select a list of the approch we can use
			% "primal"
			% "simple"
			% "hss"
			% "bb"
			% (double quotes are important) 
			for solver_approach=["bb"];%, "primal", "hss","bb"];
				disp(solver_approach)
				% for each solver approach this funciton generate a list
				% of linear solve configurations. 
				[ctrls,labels]=set_linear_algebra_ctrl(solver_approach,rec);

				
				% cycle all linear algebra controls generated with set_linear_algebra_ctrl
				for i_ctrl=1:size(ctrls,2)
					linear_algebra_ctrl=ctrls(1,i_ctrl);					
					linear_algebra_label=labels(1,i_ctrl);
					disp(linear_algebra_label)

				
					% create directories to store data
					folder=sprintf('%s/%s_rec%d/%s/',folder_runs,test_case,rec,linear_algebra_ctrl.sol);
					if ~isfolder(folder)
						mkdir(folder);
					end

					%
					% set filename for log file
					%
					if read_from_file == 0
						experiment_label = strcat(test_case_label,linear_algebra_label)
					else
						experiment_label = strcat(test_case_label,linear_algebra_label);
						pre = sprintf('restar%d_',read_from_file);
						experiment_label = strcat(pre,experiment_label)
					end
					filename = strcat(folder,experiment_label)

					% write a head with a resume of test case
					IP_ctrl.save_log = 1;
					IP_ctrl.file_log = strcat(filename,'.log');
					if exist(IP_ctrl.file_log, 'file') == 2
						delete(IP_ctrl.file_log)
					end
					logID = fopen(IP_ctrl.file_log,'w');
					fprintf(logID,'mesh type      = %d\n',mesh_type);
					fprintf(logID,'rec type       = %d\n',rec);
					fprintf(logID,'ncellphi= %d ncellrho=%d n time step=%d \n',ncell_phi,ncell_rho,Ntime);
					fclose(logID);	

					% csv file with performance resume
					% check in l2otp_solve the quantities saved
					IP_ctrl.save_csv=1;
					IP_ctrl.file_csv=strcat(filename,'.csv');

					% h5 with approximate solutions
					IP_ctrl.save_h5=0;
					IP_ctrl.file_h5=strcat(filename,'.h5');


					% Define the functional as matlab symbolic function
					% for example:
					% 
					% then use symbolic2f_df_ddf to define the triples
					% {function, derivative, second_derivate}
					% used by the l2_solve
					% `
					% syms entropy(r);
					% entropy = r*log(r);
					% f_df_ddf = symbolic2f_df_ddf(entropy)
					
					% solve with interior point
					[phi,rho,ierr,info_solver] = solve_benamoubrenier(cells,nodes,Ntime,...
																					rho_in,rho_f,...
																					1, rec, 1e-5);

					% Print errors w.r.t to exact solution
					if bc_sol==1&&compute_err==1
						% Print W2 distance
						rho_all = [rho_in;rho;rho_f];

						% rebuild some matrices
						[mass_phi,mass_edge_phi,div_phi,grad_phi] = build_matrices(grid_phi);
						gradt = assembleGradt(Ntime,ncell_phi,size(grad_phi,1),grad_phi);
						Mst = assembleMst(Ntime,size(grad_phi,1),mass_edge_phi);
						RHt = assembleRHt(Ntime,ncell_phi);
						It = assembleIt(Ntime,ncell_phi,ncell_rho,I);

						% compute W2
						%approx_W2 = compute_cost(grid_phi.ind,grid_phi.edges,grid_phi.mid,grid_phi.cc,...
						%												 gradt,Mst,RHt,It,N,rho_all,phi,rec);	
						%fprintf('%35s %1.4e \n','Approximated Wasserstein distance: ',approx_W2);
						
						
						% compare w.r.t. exact solution
						[err_cost,err_p,err_rhoth]=compute_errors(approx_W2,rho_in,rho,rho_f,phi,... % variable 
																											grid_rho,grid_phi,RHt,It,... % geometry 
																											mass,exact_W2,exact_rho_function,exact_phi_function); % exact formulas
						
						msg=sprintf('%10s %1.4e \t %11s %1.4e \t %11s %1.4e','W2-error: ',...
												err_cost,'phi-error: ',err_p,'rho-error: ',err_rhoth);
						fprintf('%s\n',msg);

						% include in log file
						logID = fopen(IP_ctrl.file_log,'a+');
						fprintf(logID,'%s\n',msg);
						fclose(logID);
					end
					
					% plot
					if (plot_figures)
						plot_rhos(grid_rho,rho_in,rho,rho_f)
					end
				end	
			end

		end  % end mesh type
	end % mesh size
end % time size
end  % test case


