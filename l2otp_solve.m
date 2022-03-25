function [phi,rho,slack,info_solver] = l2otp_solve(grid_rho, grid_phi,I, rec, Ntime,...
																									 IP_ctrl,linear_solver_ctrl,... 
																									 rho_initial,rho_final, phi_rho_slack_initial_guess)

	% SPATIAL AND TEMPORAL DISCRETIZATION
	% grid_rho (class in TPFA_grid.m) ::  grid where the problem is discreted
	% grid_phi (class in TPFA_grid.m) ::  grid where potential phi is defined
	% I (sparse matrix) ::  interpolator from grid_rho to grid_phi (it may be the identity)
	% rec (integer) :: rec==1 : linear recostruction
	%               :: rec==2 : harmonic recostruction
	% Ntime(integer):: number of time steps. Deltat=1/Ntime

	% ALGORITHM CONTROLS
	% IP_ctrl (class in IP_controls.m) :: IP algorithm controls
	% linear_solver_ctrl (structure defined in set_linear_algebra_controls.m ) :: structire
	%      containg controls for solving Newton system JF dx = - F 
	
	% PROBLEM DATA
	% rho_initial, rho_final (real arrays of grid_rho.ncell) :: non negative arrays of size grid_rho.ncell
	%                           described initial and final density distribution
	% phi_rho_slack_initial_guess (real arrays) :: initial guess of the problem

	% returns:
	% phi (real array of size grid_phi.ncell*(Ntime+1) ) with the potential 
	% rho (real array of sixe grid_rho.ncell*(Ntime)) density defined on the mesh grid_rho 
	% slack (real array of sixe grid_rho.ncell*(Ntime) ) with slackness varaible defined on the mesh grid_rho 
	% info_solve (structure created at te end) :: info solver (IP iterations, linear algebra iteration, etc)
	
	% open file to store algorithm info
	if (IP_ctrl.save_csv)
		csvID = fopen(IP_ctrl.file_csv,'w');
		fprintf(csvID,'  nrho,    np,    nt, step,    error,newton,  outer,   inner,  innernequ, minrho,    mu, theta,cpulinsys,cpuassemb,cpuprec, inner2,innernequ2,inner3,innernequ3\n');
	end


	% open log file in append mode
	if (IP_ctrl.save_log)
		logID = fopen(IP_ctrl.file_log,'a+');
	end

	% clean data 
	if (IP_ctrl.save_h5)
		if exist(IP_ctrl.file_h5, 'file')==2
			delete(IP_ctrl.file_h5);
		end
	end
	

	% Short-hand of IP parameters:
	eps_0 = IP_ctrl.eps_0; % tolerance for the final solution
	eps_mu = IP_ctrl.eps_mu; % tolerance Newton algorithm 
	k2max = IP_ctrl.k2max; % maximum number of inner (Newton) iterations
	k1max = IP_ctrl.k1max; % maximum number of outer iterations
	theta0 = IP_ctrl.theta0; % decay ratio for the perturbation parameter mu
	theta_min = IP_ctrl.theta_min;
	theta_max = IP_ctrl.theta_max;
	alpha_min = IP_ctrl.alpha_min; % minimal step length accepted
	mu0 = IP_ctrl.mu; %initial relaxation

	
	% flags for IP algorithm
	mu0_up = IP_ctrl.mu0_up;
	ierr=0;
	err_solver=0;
	err_theta=0;
	err_alpha_min=0;

	

	% dimensions
	N = Ntime;
	Nt = N + 1;
	ncell_phi = grid_phi.ncell;
	ncell_rho = grid_rho.ncell;
	nei = size(grid_phi.ind.internal,1);
	Np = (N+1)*ncell_phi; % total number of dofs for the potential
	Nr = N*ncell_rho; % total number of dofs for the density
	tne = (N+1)*nei; % total number of internal edges
	

	% set initial solution
  uk = phi_rho_slack_initial_guess; % initial condition of the barrier method
	
	% normalize rho
	for k = 1:N
		rhomass=grid_rho.area'*uk(Np+1+(k-1)*ncell_rho:Np+k*ncell_rho);
		uk(Np+1+(k-1)*ncell_rho:Np+k*ncell_rho) = uk(Np+1+(k-1)*ncell_rho:Np+k*ncell_rho)/rhomass;
	end


	% Assemble matrices associated to spatial discretization

	% local matrices
	local_ass = tic;
	[mass_phi,mass_edge_phi,div_phi,grad_phi] = build_matrices(grid_phi);
	[mass_rho,mass_edge_rho,div_rho,grad_rho] = build_matrices(grid_rho);	
	disp('local matrices assembled in')
	disp(toc(local_ass))

	% global matrices
	global_ass = tic;
	Dt = assembleDt(N,ncell_phi);
	divt = assembleDivt(N,ncell_phi,nei,div_phi);
	gradt = assembleGradt(N,ncell_phi,nei,grad_phi);
	M_phi = assembleMxt(N,ncell_phi,mass_phi);
	M_rho = assembleMxt(N,ncell_rho,mass_rho);
	Mst = assembleMst(N,nei,mass_edge_phi);
	RHt = assembleRHt(N,ncell_phi);
	It = assembleIt(N,ncell_phi,ncell_rho,I);
	I_all = [sparse(ncell_rho,Nr);speye(Nr,Nr);sparse(ncell_rho,Nr)];

	if rec==1
    Rs=Ktos2D(grid_phi.ind,grid_phi.edges,...
							grid_phi.cc,grid_phi.mid);
    Rst = repmat({Rs},1,N+1);
    Rst = blkdiag(Rst{:});
    %clear Rs
		Rs_rho  = Ktos2D(grid_rho.ind,grid_rho.edges,grid_rho.cc,grid_rho.mid);
    Rst_rho = repmat({Rs_rho},1,N);
    Rst_rho = blkdiag(Rst_rho{:});
	else
    Rst = [];
	end

	nei_rho = size(div_rho,2);
	divt_rho  = assembleDivt(N-1,grid_rho.ncell,nei_rho,div_rho); 
	gradt_rho = assembleGradt(N-1,grid_rho.ncell,nei_rho,grad_rho);

	
	
	disp('global matrices assembled in ')
	disp(toc(global_ass))

	% Preprocess finished
	
	
	
	
	% Start IP algorithm
	itk1 = 0; % counter of IP iterations
	tit = 0; % counter of the total number of Newton iterations
	theta = theta0;
	mu = mu0/theta;
	delta_0 = 2*eps_0;

	% cpu-time record
	cpu_total=0.0;
	cpu_assembly=0.0;
	cpu_linsys=0.0;

	% IP cycle 
	while true
		assembly=tic;
		total=tic;
		if delta_0 < eps_0
			break
		end
		itk1 = itk1+1;
		if itk1>k1max
			% set IP controls
			if mu0_up==1
				mu0 = 5*mu0;
				theta = theta0;
				s0 = mu0./rho0;
				mu = mu0/theta;
				itk1 = 0;
				tit = 0;
				uk = [phi0;rho0;s0];
				delta_0 = 2*eps_0;
				continue
			else
				break
			end
		end
    
    mu = theta*mu;
    %eps_mu = eps_0;
    itk2 = 0;
    flag2 = 0;

		resvar = set_resval;
		
		% save before update
    if (IP_ctrl.save_h5 > 0)
      data_name=sprintf('/DS%d',itk1);
      h5create(IP_ctrl.file_h5,data_name,[Np+2*Nr+2])
      h5write(IP_ctrl.file_h5,data_name,[uk;mu;theta]')
    end
		% save before next Newton cycle
		uk_before=uk;


    % reset counter linear algebra
    sum_assembly=0;
    sum_total=0;
    sum_linsys=0;
    sum_prec=0;
    sum_iter_newton=0;
    sum_iter_outer_linear=0;
    sum_iter_inner_linear=0;
    sum_iter2_inner_linear=0;
    sum_iter3_inner_linear=0;

		% start Newton cycle
    while true
      total=tic;
      

      % assembly non-linear equations
			assembly=tic;
      [rhosk]=compute_rhosigma(grid_phi.ind,grid_phi.edges,grid_phi.cc,grid_phi.mid,...
															 N,rho_final,rho_initial,gradt,Mst,RHt,It,Rst,rec,uk,'rhos');
      [drhosk]=compute_rhosigma(grid_phi.ind,grid_phi.edges,grid_phi.cc,grid_phi.mid,...
																N,rho_final,rho_initial,gradt,Mst,RHt,It,Rst,rec,uk,'drhos');

      ctime=tic;
			% note the we relax the rho_in and rho_f. 
      OC = Fkgeod(N,(rho_final+mu)/(1+mu),(rho_initial+mu)/(1+mu),...
									Dt,divt,M_phi,M_rho,Mst,gradt,It,rhosk,drhosk,uk,mu);
      FOCtime=toc(ctime);
      delta_mu = norm([OC.p;OC.r;OC.s]);

			% print info on optimality condition
      state_message=sprintf('%d - |OC.p|=%1.4e |OC.r|=%1.4e |OC.s|=%1.4e - CPU %1.4e' ,...
														itk2+1, norm(OC.p),norm(OC.r),norm(OC.s),FOCtime);		
      if (IP_ctrl.verbose >= 2)
				fprintf('%s \n',state_message);
      end
			if (IP_ctrl.save_log)
				fprintf(logID,'%s \n',state_message);
			end
			

			% check convergence or adjust IP adjust IP parameters 
      if delta_mu < eps_mu	
        if itk2<4
          if 0.8*theta>=theta_min
            theta = 0.8*theta;
          else 
            theta = theta_min;
          end
        end
        tit = tit+itk2;
        break
      end
      itk2 = itk2+1;
      if itk2 > k2max || flag2==1
				% if the step length is too small, or the number of Newton iterations 
				% exceeds the maximum, repeat the iteration with a bigger mu
        if itk1==1
					ierr = 1;
          itk1 = k1max;
          break
        end
        if theta+0.2*(1-theta)<=theta_max
          mu = mu/theta;
          theta = theta+0.2*(1-theta);
          mu = theta*mu;
        else
					ierr = 2;
          break
        end
        uk = [phimu;rhomu;smu];
        tit = tit+itk2;
        itk2 = 0;
        flag2 = 0;
        continue
      end
			cpu_assembly=toc(assembly);
      
      
      
      
		  % Compute the jacobian of the system of equations
			ctime=tic;
      [ddrhosak]=compute_rhosigma(grid_phi.ind,grid_phi.edges,grid_phi.cc,grid_phi.mid,...
																	N,rho_final,rho_initial,gradt,Mst,RHt,It,Rst,rec,uk,'ddrhosa');
      JOC = JFkgeod(N,Dt,divt,M_phi,M_rho,gradt,It,rhosk,drhosk,ddrhosak,...
										uk,I,Rs,Rst_rho,divt_rho,gradt_rho);

      sum_assembly=sum_assembly+toc(assembly);
      JFOCtime=toc(ctime);
      if ( IP_ctrl.verbose >=3 )
        fprintf('CPU ASSEMBLY: TOTAL %1.4e - FOC=%1.4e -JOC=%1.4e \n',toc(assembly),FOCtime,JFOCtime)
      end

			resvar.set(IP_ctrl.kel,IP_ctrl.min_outer_tol,delta_mu,Np+2*Nr);
			ctrl_outer.tolerance=resvar.etak;

			% solve linear system J d = -F 
			timelinsys=tic;
			[omegak, linear_solver_info,norm_ddd,resume_msg] = solvesys(JOC,OC, linear_solver_ctrl);
			
			% print linear solver info
			if (strcmp(resume_msg,''))
				print_info_solver(linear_solver_info)
				print_info_solver(linear_solver_info,logID)
			else
				fprintf('%s\n',resume_msg);
				fprintf(logID,'%s\n',resume_msg);
			end
			
			if (~ (linear_solver_info.flag == 0) )
				disp('ERROR')
				fprintf(logID,'%s\n','ERROR');
				if (linear_solver_info.relres > 5*resvar.etamax)
					ierr=3;
					uk=uk_before;
					break
				end
			end

			%
			% collect data ad save it to csv file
			%

			% measure linear solver performance		
			sum_linsys = sum_linsys+toc(timelinsys); % cumulative cpu time in iterative solver
			sum_total  = sum_total+toc(total);  
			sum_prec   = sum_prec+linear_solver_info.prec_cpu; % cumulative cpu-time in preprocess before linear solver

			sum_iter_outer_linear  = sum_iter_outer_linear + linear_solver_info.outer_iter; 
			sum_iter_inner_linear  = sum_iter_inner_linear + linear_solver_info.inner_iter;
			sum_iter2_inner_linear = sum_iter2_inner_linear + linear_solver_info.inner_iter2;
			sum_iter3_inner_linear = sum_iter3_inner_linear + linear_solver_info.inner_iter3;

			cpu_total=cpu_total+sum_total;
			cpu_assembly=cpu_assembly+sum_assembly;
			cpu_linsys=cpu_linsys+sum_linsys;
			
			
			
			% Linesearch just to ensure that rho and s stay positive
			alpha_k = 1;
			while any(uk(Np+1:end)+alpha_k*omegak(Np+1:end)<=0)
				alpha_k = 0.5*alpha_k;
				if alpha_k < alpha_min
					alpha_k = alpha_min;
					ierr=4;
					uk=uk_before;
					break
				end
			end

			uk = uk + alpha_k*omegak;
			
			if alpha_k==alpha_min
				disp('Step too small')
				flag2 = 1;
			end

			% se phi to sum to zero
			%uk(1:Np)=ort_proj(uk(1:Np),ones(1,Np));
			
			% normalize rho
			for k = 1:N
				rhomass=grid_rho.area'*uk(Np+1+(k-1)*ncell_rho:Np+k*ncell_rho);
				%fprintf('mass rho(%d)-1.0=%1.1e\n',k,rhomass-1.0)
				uk(Np+1+(k-1)*ncell_rho:Np+k*ncell_rho) = uk(Np+1+(k-1)*ncell_rho:Np+k*ncell_rho)/rhomass;
			end

			if (IP_ctrl.save_h5 > 1 )
				data_name=sprintf('/DS%d_inner%d',itk1,itk2);
				h5create(IP_ctrl.file_h5,data_name,[Np+2*Nr+2])
				h5write(IP_ctrl.file_h5,data_name,[uk;mu;theta]')
			end

			
		end

    % exit from IP cy
    if ierr ~= 0;
			uk=uk_before;
      break
    end
    
    phimu = uk(1:Np);
    rhomu = uk(Np+1:Np+Nr);
    smu = uk(Np+Nr+1:end);

    
    
    % error bound on optimality
    delta_0 = (N/Nt)*sum(grid_phi.area)*mu;

    % message resuming IP solver status
    state_message=sprintf('%5s %3i %4s %1.4e %7s %3i  %10s %1.4e %4s %1.2e %7s %1.2e',...
													'IP it:',itk1,...
													' Err',delta_0,...
													' Newton it.', itk2,...
	    										" min(rho) ",min(uk(Np+1:Np+Nr))," mu ",mu," theta ",theta);
    if ( IP_ctrl.verbose>0 )
      fprintf('%s \n', state_message)	
		end

		if (IP_ctrl.verbose >=1 )
			fprintf(logID,'%s \n',state_message);
		end
		
	
	
		
		if (IP_ctrl.save_csv)
			fprintf(csvID,'%6d,%6d,%6d,%5d,%8.3e,%6d,%7d,%8d,%11d,%1.1e,%1.1e,%1.1e,%1.1e,%1.1e,%1.1e,%8d,%11d,%8d,%11d\n',...
							ncell_rho,ncell_phi,N,...
							itk1,delta_0,itk2,sum_iter_outer_linear,uint64(sum_iter_inner_linear),...
							linear_solver_info.inner_nequ,min(uk(Np+1:Np+Nr)),mu,theta,...
							sum_linsys,sum_assembly,sum_prec,...
							uint64(sum_iter2_inner_linear),linear_solver_info.inner_nequ2,...
							uint64(sum_iter3_inner_linear),linear_solver_info.inner_nequ3);
		end
	end

	% save last last
	if ierr== 0 && IP_ctrl.save_h5 > 0
		data_name=sprintf('/DS%d',itk1+1);
		h5create(IP_ctrl.file_h5,data_name,[Np+2*Nr+2])
		h5write(IP_ctrl.file_h5,data_name,[uk;mu;theta]')
	end

	% print resume messagge
	fprintf('%17s %4i %7s %1.4e \n','Total Newton its:',tit,'Error: ',delta_0)
	fprintf(logID,'%17s %4i %7s %1.4e \n','Total Newton its:',tit,'Error: ',delta_0);
	fprintf(logID,'%19s %1.4e %21s %1.4e \n','Total linsys time: ',cpu_linsys,'Total assembly time: ',cpu_assembly);

	
	% close files
	fclose(logID);
	fclose(csvID);


	% get the solution 
	phi   = uk(1:Np);
	rho   = uk(Np+1:Np+Nr);
	slack = uk(Np+Nr+1:Np+2*Nr);

	% set info_solver
	info_solver=struct('ierr',ierr,...
										 'ip_iter',itk1,...
										 'newton_iter',tit);

end
