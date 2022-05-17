clear all
close all


% set test case fixing intitial and final density.
% See function boundary_bc for the available options.
test_cases = ["sin","gauss","compression"];
for  test_case = test_cases;
	test_case = test_case;
	disp(test_case)

folder_runs='FTNB_runs';
if ~isfolder(folder_runs)
	mkdir(folder_runs);
end

folder_restart='initial_solution';

%
% set some globals controls for reproducing experiments
%
read_from_file = 0;
plot_figures = 0;
compute_err = 1;


%
% SPATIAL DISCRETIZATION
%
for mesh_type = 5;
	for h_i = 1:4;
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
		[grid_rho, grid_phi, I] = init_grids(mesh_type, h_i);

		%
		% TEMPORAL DISCRETIZATION delta=1/N
		%
		for dt_i = 2:5
			N=4*(2^(dt_i-1));

			% set problem dimension
			ncell_phi=grid_phi.ncell;
			ncell_rho=grid_rho.ncell;
			Np=(N+1)*ncell_phi;
			Nr=(N)*ncell_rho;


			%
			% PROBLEM DATA
			%
			% Set initial and final rho density
			% In same cases the exact solution is known
			[rho_in,rho_f,...
			 mass,midpoint,exact_phi,exact_rho,exact_W2,bc_sol] = ...
			bc_density(test_case,grid_rho.cc,grid_rho.area);

			


			test_case_label = sprintf('%s_rec%d_mesh%d_h%d_N%0.5d_',test_case,rec,mesh_type,h_i,N);


						%
			% ALGORITHM CONTROLS
			%
			% Init Interior Point controls.
			% Defalt values are set in IP_controls class
			IP_ctrl=IP_controls;
			
			% Set initial data
			% Default (phi, rho, s ) = 1,1,mu0
			if ( read_from_file == 0)
				mu0 = 1;
				phi0 = ones(Np,1);
				rho0 = (mass/sum(grid_rho.area))*ones(Nr,1);
				s0 = mu0./rho0;
				uk = [phi0;rho0;s0]; 
			else
				% this option can be used to restart the IP algorithm form a certain solution
				% str_folder=sprintf('%s_mesh%d_rec%d',test_case,mesh_type,rec)
				folder=sprintf('%s/%s_rec%d/initial_solution/',folder_runs,test_case,rec);
				str_test=sprintf('%s_rec%d_mesh%d_h%d_N%0.5d_',test_case,rec,mesh_type,h_i,N);
				h5_file2read = strcat(folder,...
															str_test(1:strlength(str_test)-1),'.h5')
				data_name=sprintf('/DS%d',read_from_file);
				data = h5read(h5_file2read,data_name);
				phi0 = data(1:Np);
				rho0 = data(1+Np:Np+Nr);
				s0   = data(1+Np+Nr:Np+2*Nr);
				mu0  = data(Np+2*Nr+1);
				theta0 = data(Np+2*Nr+2);
				%fprintf('norm initial %1.2e\n',norm(get_slice(phi0,1,ncell_phi,N+1)))
				uk = [phi0;rho0;s0];
				test_case_label=sprintf('restart%d_%s', read_from_file,test_case_label)
				IP_ctrl.mu=mu0;
				IP_ctrl.theta0=theta0;
			end

			


			for solver_approach=[11];
				% set controls
				[ctrls,labels]=set_linear_algebra_ctrl(solver_approach,rec);

				%labels

				
				% cycle all linear algrabra controls genreted with set_linear_algebra_ctrl
				%size(ctrls,2)
				
				for i_ctrl=1:size(ctrls,2)
					linear_algebra_ctrl=ctrls(1,i_ctrl);
					
					linear_algebra_label=labels(1,i_ctrl);
					%disp(test_case_label)
					disp(linear_algebra_label)

				
					
					% create directories to store data
					folder=sprintf('%s/%s_rec%d/sol%d/',folder_runs,test_case,rec,linear_algebra_ctrl.sol);
					if ~isfolder(folder)
						mkdir(folder);
					end

					%
					% set filename for files with discrtetizatio and liear algebra info
					%
					if read_from_file == 0
						experiment_label = strcat(test_case_label,linear_algebra_label)
					else
						experiment_label = strcat(test_case_label,linear_algebra_label);
						pre = sprintf('restar%d_',read_from_file);
						experiment_label = strcat(pre,experiment_label)
					end
					filename = strcat(folder,experiment_label)

					
					
					% log file, writing a small resume of test case
					IP_ctrl.save_log = 1;
					IP_ctrl.file_log = strcat(filename,'.log');
					if exist(IP_ctrl.file_log, 'file') == 2
						delete(IP_ctrl.file_log)
					end
					logID = fopen(IP_ctrl.file_log,'w');
					fprintf(logID,'mesh type      = %d\n',mesh_type);
					fprintf(logID,'rec type       = %d\n',rec);
					fprintf(logID,'ncellphi= %d ncellrho=%d n time step=%d \n',ncell_phi,ncell_rho,N);
					fclose(logID);	

					% csv file with performance resume
					IP_ctrl.save_csv=1;
					IP_ctrl.file_csv=strcat(filename,'.csv');

					% h5 with approximate solutions
					IP_ctrl.save_h5=0;
					IP_ctrl.file_h5=strcat(filename,'.h5');

					
					% solve with interior point
					[phi,rho,slack,info_solver] = ...                   % solver output and info
					l2otp_solve(grid_rho, grid_phi,I, rec, N,...        % discretization 
											IP_ctrl,linear_algebra_ctrl,...         % controls
											rho_in,rho_f, uk(1:Np+Nr+Nr)); % inputs
					
					
					% Print errors w.r.t to exact solution
					if bc_sol==1&&compute_err==1
						% Print W2 distance
						rho_all = [rho_in;rho;rho_f];

						% rebuild some matrices
						[mass_phi,mass_edge_phi,div_phi,grad_phi] = build_matrices(grid_phi);
						gradt = assembleGradt(N,ncell_phi,size(grad_phi,1),grad_phi);
						Mst = assembleMst(N,size(grad_phi,1),mass_edge_phi);
						RHt = assembleRHt(N,ncell_phi);
						It = assembleIt(N,ncell_phi,ncell_rho,I);

						% compute W2
						approx_W2 = compute_cost(grid_phi.ind,grid_phi.edges,grid_phi.mid,grid_phi.cc,...
																		 gradt,Mst,RHt,It,N,rho_all,phi,rec);	
						fprintf('%35s %1.4e \n','Approximated Wasserstein distance: ',approx_W2);
						
						
						% compare w.r.t. exact solution
						[err_cost,err_p,err_rhoth]=compute_errors(approx_W2,rho_in,rho,rho_f,phi,... % variable 
																											grid_rho,grid_phi,RHt,It,... % geometry 
																											mass,exact_W2,exact_rho,exact_phi); % exact formulas
						
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


