%include agmg folder
addpath(genpath('./agmg/'));

%% Code starts
if mesh_type == 1
    mesh_name = strcat('meshes/tri2_mesh',num2str(h_i));
elseif mesh_type == 2
    mesh_name = strcat('meshes/subtri2_mesh',num2str(h_i));
elseif mesh_type == 3
  mesh_name = strcat('meshes/sq_mesh',num2str(h_i));
elseif mesh_type == 4
      mesh_name = strcat('meshes/tri2ord_mesh',num2str(h_i));
elseif mesh_type == 5
  mesh_name = strcat('meshes/subtri2ord_mesh',num2str(h_i));
elseif mesh_type == 6
  mesh_name = strcat('meshes/sqord_mesh',num2str(h_i));
else
  disp('Mesh not supprted')
  return
end

cpu_total=0.0;
cpu_assembly=0.0;
cpu_linsys=0.0;

% load mesh
clear 'I'
load(mesh_name)


if exist('I','var')==0
    % if the injection operator I does not exist in the mesh structure,
    % the nested sub-mesh coincides with the outer one (labeled with 2h)
    ncell2h = ncell;
    nodes2h = nodes;
    cells2h = cells;
    edges2h = edges;
    ind2h = ind;
    area2h = area;
    cc2h = cc;
    mid2h = mid;
    h2h = h;
    I = speye(ncell,ncell);
end

if (~exist('approach_string','var') )
  fprintf('Set varible: ( approach_string ) with linear solver approach resume\n')
  return
end

%if (~exist('controls_string','var') )
%  fprintf('Set varible: ( controls_string ) with linear solver approach resume\n')
%  return
%end

controls.swap_sign=1;


controls.indc=(N+1)*ncell;


solver_approach=controls.sol
folder_approach=sprintf('sol%d/',solver_approach);
filename=strcat('runs/',folder_approach,str_test,approach_string);%,controls_string);
disp(filename)
controls.basename=filename;
log_filename=strcat(filename,'.log')
logID = fopen(log_filename,'w');
controls.logID=logID;

eigen_fn=sprintf('eigenS%d.csv',read_from_file);
eigen_fid = fopen(eigen_fn, 'a+');

fprintf(logID,'mesh type       = %d\n',mesh_type);
fprintf(logID,'rec type       = %d\n',rec);
fprintf(logID,'ncellphi= %d ncellrho=%d n time step=%d \n',ncell,ncell2h,N);
fprintf(logID,'solver approach = %d\n',solver_approach);
fprintf(logID,'string approach = %s\n',approach_string);
fprintf(logID,'outer solver :\n');
ctrl_outer.info(logID);
fprintf(logID,'solver block 11:\n');
ctrl_inner11.info(logID);
fprintf(logID,'solver block 22: \n');
ctrl_inner22.info(logID);



csv_filename=strcat(filename,'.csv')
csvID = fopen(csv_filename,'w');
fprintf(csvID,'  nrho,    np,    nt, step,    error,newton,  outer,   inner,  innernequ,  minrho,      mu,   theta, cpulinsys, cpuassemb,   cpuprec, inner2,innernequ2,inner3,innernequ3\n');

% Boundary conditions:
[rho_in,rho_f,mass,midpoint,potential,geodesic,W2,bc_sol] = bc_density(test_case,cc2h,area2h);

% % Barrier method's parameters:
% eps_0 = 1e-6; % tolerance for the final solution
% k2max = 30; % maximum number of inner (Newton) iterations
% k1max = 20; % maximum number of outer iterations
% theta0 = 0.2; % decay ratio for the perturbation parameter mu
% theta_min = 0.2;
% theta_max = 0.2;
% alfamin = 0.1; % minimal step length accepted
% mu0_up = 0;
flag_theta = 0;


tnp = (N+1)*ncell; % total number of dofs for the potential
tnr2h = N*ncell2h; % total number of dofs for the density
tne = (N+1)*nei; % total number of internal edges


% set initial data
if ( read_from_file>0)
  data_name=sprintf('/DS%d',read_from_file);
  data = h5read(h5_file2read,data_name);
  phi0 = data(1:tnp);
  rho0 = data(1+tnp:tnp+tnr2h);
  s0   = data(1+tnp+tnr2h:tnp+2*tnr2h);
  mu0  = data(tnp+2*tnr2h+1);
  theta0 = data(tnp+2*tnr2h+2);
  fprintf('norm initial %1.2e\n',norm(get_slice(phi0,1,ncell,N+1)))
  uk = [phi0;rho0;s0];
else
  mu0 = 1;
  phi0 = ones(tnp,1);%zeros(tnp,1);
  rho0 = (mass/sum(area))*ones(tnr2h,1);
  s0 = mu0./rho0;
  uk = [phi0;rho0;s0]; % initial condition of the barrier method
end

Np = ncell*Nt;
Nr = ncell2h*N;

% normalize rho
for k = 1:N
	rhomass=area2h'*uk(Np+1+(k-1)*ncell2h:Np+k*ncell2h);
  uk(Np+1+(k-1)*ncell2h:Np+k*ncell2h) = uk(Np+1+(k-1)*ncell2h:Np+k*ncell2h)/rhomass;
end

save_data=controls.save_data;
if (save_data)
  filename_h5=([strcat('runs/',folder_approach,'PhiRhoSMuTheta',str_test,approach_string,'.h5')]);
  if exist(filename_h5, 'file')==2
    delete(filename_h5);
  end
end


% Assemble matrices
local_ass = tic;
Mx = spdiags(area,0,ncell,ncell);
Mx2h = spdiags(area2h,0,ncell2h,ncell2h);
ds = edges(ind.internal,5).*edges(ind.internal,6);
Ms = spdiags(ds,0,nei,nei);
div = Div2D(ncell,nei,ind,edges); % divergence matrix
grad = -Ms\div'; % gradient matrix

nei2h=size(ind2h.internal,1)
ds_rho = edges2h(ind2h.internal,5).*edges2h(ind2h.internal,6);
size(ds_rho)
Ms_rho = spdiags(ds_rho,0,nei2h,nei2h);
div_rho = Div2D(ncell2h,nei2h,ind2h,edges2h); % divergence matrix
grad_rho = -Ms_rho\div_rho'; % gradient matrix
divt_rho  = assembleDivt(N-1,ncell2h,nei2h,div_rho); 
gradt_rho = assembleGradt(N-1,ncell2h,nei2h,grad_rho);




disp('local matrices assembled in')
disp(toc(local_ass))

% global matrices
global_ass = tic;
Dt = assembleDt(N,ncell);
divt = assembleDivt(N,ncell,nei,div);
gradt = assembleGradt(N,ncell,nei,grad);
Mxt = assembleMxt(N,ncell,Mx);
Mxt2h = assembleMxt(N,ncell2h,Mx2h);
Mst = assembleMst(N,nei,Ms);
RHt = assembleRHt(N,ncell);
It = assembleIt(N,ncell,ncell2h,I);
I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
if rec==1
    Rs=Ktos2D(ind,edges,cc,mid);
    Rst = repmat({Rs},1,N+1);
    Rst = blkdiag(Rst{:});
    %clear Rs
		Rs_rho  = Ktos2D(ind2h,edges2h,cc2h,mid2h);
    Rst_rho = repmat({Rs_rho},1,N);
    Rst_rho = blkdiag(Rst_rho{:});
else
    Rst = [];
end
disp('global matrices assembled in ')
disp(toc(global_ass))

itk1 = 0;
tit = 0; % counter of the total number of Newton iterations
theta = theta0;
mu = mu0/theta;
delta_0 = 2*eps_0;

masses=zeros(N,1);
masses_increment=zeros(N,1);
imbalance_res_phi=zeros(Nt,1);
imbalance_increment_phi=zeros(Nt,1);

Np = ncell*Nt;
Nr = ncell2h*N;






while true
  assembly=tic;
  total=tic;
  if delta_0 < eps_0
    break
  end
  itk1 = itk1+1;
  if itk1>k1max
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
      disp(mu0_up)
      break
    end
  end
    
    mu = theta*mu;
    eps_mu = eps_0;
    itk2 = 0;
    flag2 = 0;

    resvar = set_resval;
    
    if (save_data)
      data_name=sprintf('/DS%d',itk1);
      h5create(filename_h5,data_name,[Np+2*Nr+2])
      h5write(filename_h5,data_name,[uk;mu;theta]')
    end



    
    sum_assembly=0;
    sum_total=0;
    sum_linsys=0;
    sum_prec=0;
    sum_iter_newton=0;
    sum_iter_outer_linear=0;
    sum_iter_inner_linear=0;
    sum_iter2_inner_linear=0;
    sum_iter3_inner_linear=0;
    
    while true
      total=tic;
      assembly=tic;
      
      % message="BEGIN ASSEMBLY FOC";
      % if verb>1
      % 	if (  itk2 == 0)
      % 	  state_message=sprintf(' \n')	;		
      % 	end	
      % 	fprintf('%s \n',message);
      % end

      [rhosk]=compute_rhosigma(ind,edges,cc,mid,N,rho_f,rho_in,gradt,Mst,RHt,It,Rst,rec,uk,'rhos');
      [drhosk]=compute_rhosigma(ind,edges,cc,mid,N,rho_f,rho_in,gradt,Mst,RHt,It,Rst,rec,uk,'drhos');

      ctime=tic;
      OC = Fkgeod(N,(rho_f+mu)/(1+mu),(rho_in+mu)/(1+mu),Dt,divt,Mxt,Mxt2h,Mst,gradt,It,rhosk,drhosk,uk,mu);
      FOCtime=toc(ctime);
      delta_mu = norm([OC.p;OC.r;OC.s]);
      
      state_message=sprintf('%d - |OC.p|=%1.4e |OC.r|=%1.4e |OC.s|=%1.4e - CPU %1.4e' ,...
			    itk2+1, norm(OC.p),norm(OC.r),norm(OC.s),FOCtime);
      fprintf(logID,'%s \n',state_message);
      if verb>1
	fprintf('%s \n',state_message);
      end
      %fprintf(2,'%s \n',state_message);


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
          itk1 = k1max;
          break
        end
        if theta+0.2*(1-theta)<=theta_max
          mu = mu/theta;
          theta = theta+0.2*(1-theta);
          mu = theta*mu;
        else
          flag_theta = 1;
          break
        end
        uk = [phimu;rhomu;smu];
        tit = tit+itk2;
        itk2 = 0;
        flag2 = 0;
        continue
      end
      
      sys_name=sprintf('out%02d_in%02d',itk1+1,itk2+1);
      controls.sys_name=sys_name;
      
		     % Compute the jacobian of the system of equations
      % message="BEGIN ASSEMBLY JFOC";
      % if verb>1
      %   fprintf('%s \n',message)
      % end
      ctime=tic;
      % Compute the jacobian of the system of equations
      [ddrhosak]=compute_rhosigma(ind,edges,cc,mid,N,rho_f,rho_in,gradt,Mst,RHt,It,Rst,rec,uk,'ddrhosa');
      
      JOC = JFkgeod(N,Dt,divt,Mxt,Mxt2h,gradt,It,rhosk,drhosk,ddrhosak,...
										uk,I,Rs,Rst_rho,divt_rho,gradt_rho);

      if (augmented)
        % aug. Jacobian need the original OC
				[OC,JOC] = newton_augmentation(uk,OC,JOC,option,line);
      end
      sum_assembly=sum_assembly+toc(assembly);
      JFOCtime=toc(ctime);
			
      if verb>2
        fprintf('CPU ASSEMBLY: TOTAL %1.4e - FOC=%1.4e -JOC=%1.4e \n',toc(assembly),FOCtime,JFOCtime)
      end

			
			if ( compute_eigen)
				[A,B1T,B2,C,f1,f2,W_mat]=get_saddle_point(JOC,OC);
				B1T_time=-JOC.B1T_time;
				B1T_space=-JOC.B1T_space;

				P = @(y) P_apply(y,area2h);
				PT = @(y) PT_apply(y,area2h);
				
				if (times_H)
					% A,   B1T*mat_H';
  				% mat_H*B2, -mat_H*C*mat_H'; 
					[A,B1T_time,B1T_space,B2,C,f1,f2,H_mat] = times_H(A,B1T_time,B1T_space,B2,C,f1,f2,JOC)
				end

				if (ground)
					% ground system 
					[A,B1T_time,B1_space,B2_time,B2_space,f1]=ground_saddle(A,B1T_time,B1T_space,f1,N);
					B1T=B1T_time+B1_space;
					B2 =B2_time +B2_space;
				end

				
				if (diagonal_scaling)
					% scale system by its diagonal
					[sqrt_diagA,sqrt_diagC,inv_sqrt_diagA,inv_sqrt_diagC] = get_diagonal_scaling (A,C);
					
					[A,B1T,B2,C,f1,f2,W_mat]=scale_saddle_point(A,B1T,B2,C,...
																											f1,f2,W_mat,...
																											inv_sqrt_diagA,...
																											inv_sqrt_diagC);
					B1T_time=inv_sqrt_diagA* B1T_time *inv_sqrt_diagC;
					B1T_spce=inv_sqrt_diagA* B1T_space *inv_sqrt_diagC;
				end

				build_real_schur=1;
				if (build_real_schur)
					if (ground)
						[Stt,Stx,Sxx]=assemble_dual_schur_components(A,size(B1T_time,2),N,...
																												 @(y) B1T_time*(y),...
																												 @(y) B1T_space*(y),...
																												 @(x) P(B2_time*x),...
																												 @(x) P(B2_space*x));
						S=Stt+Stx+Sxx+Stx';
					else
						if (times_H)
							[Stt,Stx,Sxx]=assemble_dual_schur_components(A,size(B1T_time,2),N,...
																													 @(y) B1T_time*y,...
																													 @(y) B1T_space*y,...
																													 @(x) B1T_time'*x,...
																													 @(x) B1T_space'*x);
							
						else
							[Stt,Stx,Sxx]=assemble_dual_schur_components(A,size(B1T_time,2),N,...
																													 @(y) B1T_time*PT(y),...
																													 @(y) B1T_space*PT(y),...
																													 @(x) B1T_time'*x,...
																													 @(x) B1T_space'*x);
							
							CPT=matrix_times_prec(C,PT);
							PCPT=apply_prec_matrix(P,CPT);
							Stt=prec_times_matrix(P,Stt);
							Stx=prec_times_matrix(P,Stx);
							Sxx=prec_times_matrix(P,Sxx);
						end
						Schur=Stt+Sxx+Stx+Stx';
					end
				end

				%
				if (0)
					figure
					[real,imag,eigenvec,perm]= study_eigenvalues(Schur,'never gonna works');
					plot(real,imag,'o')
					
					figure
					semilogy(real)
					title('S')

					inverse_S = sparse_inverse;
					ctrl_loc = ctrl_solver;
					index_agmg=1;
					ctrl_loc.init('direct',...
    										1e-13,...
    										ctrl_inner11.itermax,...
    										ctrl_inner11.omega,...
    										1,...
    										'S',...
												index_agmg);
					
					inverse_S.init(Schur+1e-12*speye(Nr,Nr),ctrl_loc);
					inv_S=@(y) inverse_S.apply(y);

					C_S=prec_times_matrix(inv_S,PCPT);
					Stt_S=prec_times_matrix(inv_S,Stt);
					Sxx_S=prec_times_matrix(inv_S,Sxx);
					Smixed_S=prec_times_matrix(inv_S,Stx+Stx');

					[real,imag,eigenvec,perm]= study_eigenvalues(C_S,'never gonna works');
					figure
					plot(real,imag,'o')
					title('C_S')

					[real,imag,eigenvec,perm]= study_eigenvalues(Stt_S,'never gonna works');
					figure
					plot(real,imag,'o')
					title('Stt_S')

					
					[real,imag,eigenvec,perm]= study_eigenvalues(Sxx_S,'never gonna works');
					figure
					plot(real,imag,'o')
					title('Sxx_S')

					[real,imag,eigenvec,perm]= study_eigenvalues(Smixed_S,'never gonna works');
					figure
					plot(real,imag,'o')
					title('Sxx_S')
								
					return
				end
				
				if (1)
					% Pspace=assemble_space_projector(N+1,I');  
					% Ptime=assemble_time_projector(N,ncell2h);
					% pr = Pspace'*Ptime';
					
					% %time_Laplacian=(Nt*Dt*It*I_all)'*Mxt*Mxt*(Nt*Dt*It*I_all);
					inv_Mxt  = sparse(1:Np,1:Np,(1./spdiags(Mxt,0))',Np,Np);
					% inv_Mrho = sparse(1:Nr,1:Nr,(-1./spdiags(JOC.rs))',Nr,Nr);
					
					% A_rho = pr'*inv_Mxt*A*inv_Mxt* pr;
					% A_rho = Ptime*Pspace*inv_Mxt*A*inv_Mxt* Pspace'*Ptime'*inv(Ptime*Ptime');

					% %
					% %	(C+B2 invA B1^T)^{-1}=
					% % (C+B2*B1T*invA)^{-1}=
					% % (A_rho)*(C*A_rho+time_laplacian)^{-1}=:prec
					% % 
					% approx_S=C*A_rho+B2*inv_Mxt*inv_Mxt*B1T;%time_Laplacian;
					lrb='b'
					[pr,A_rho,approx_S,start]=build_approx_Schur_components(A,B1T,B2,C,lrb,JOC);
					
					inverse_SC=sparse_inverse;
					ctrl_loc=ctrl_solver;
					ctrl_loc.init('agmg',...
					 							1e-1,...
					 							200,...
					 							0,...
					 							1,...
					 							'invS_approx',...
					 							84);
					inverse_SC.init(approx_S+1e-12*speye(Nr),ctrl_loc);
				
					inverse_approx_S = @(y) ...
															compute_inverse_approx_S(y,...
																											 @(z) inverse_SC.apply(z),...
																											 A_rho,pr,...
																											 lrb,start,Nr,size(approx_S,1))

					%inverse_approx_S = @(y) inverse_SC.apply(y)
					
				  % Define action of schur complement
					%
					inverseA=sparse_inverse;
					ctrl_loc=ctrl_solver;
					index_agmg=1;
					ctrl_loc.init('direct',...
    										1e-13,...
    										ctrl_inner11.itermax,...
    										ctrl_inner11.omega,...
    										0,...
    										'A',...
												index_agmg);
					
					inverseA.init(A+1e-13*speye(Np,Np),ctrl_loc);
					inverseA.cumulative_iter=0;
					inverseA.cumulative_cpu=0;

					% define function
					invA = @(x) inverseA.apply(x);

					if(0)
						S_tilde=zeros(Nr,Nr);
						temp=zeros(Nr,1);
						inv_Mxt  = sparse(1:Np,1:Np,(1./spdiags(Mxt,0))',Np,Np);

						Q1=inv_Mxt;
						Q2=inv_Mxt;
						
						P = @(y) P_apply(y,area2h);
						PT = @(y) PT_apply(y,area2h);
						BBT=prec_times_matrix(P,matrix_times_prec(B2*Q1*B1T,PT));

						inverse_BBT=sparse_inverse;
						ctrl_loc=ctrl_solver;
						index_agmg=1;
						ctrl_loc.init('direct',...
    											1e-13,...
    											ctrl_inner11.itermax,...
    											ctrl_inner11.omega,...
    											1,...
    											'BBT',...
													index_agmg);
						inverse_BBT.init(BBT+1e-13*speye(Nr,Nr),ctrl_loc);

						lsc= @(y) inverse_BBT.apply(P(B2*Q1*A*Q2*(B1T*PT(inverse_BBT.apply(y)))));

						%S=P(Stt+Sxx+Stx+Stx');
						temp=prec_times_matrix(lsc,Stt+Stx+Stx'+Sxx);

						% %ragmgeturnagmg
						%figure
						[real,imag,eigenvec,perm]= study_eigenvalues(temp,'never gonna works');
						figure
						plot(real,imag,'o')

						figure
						semilogy(real)
						return

					end

					if(0)
						S_tilde=zeros(Nr,Nr);
						temp=zeros(Nr,1);
						inv_Mxt  = sparse(1:Np,1:Np,(1./spdiags(Mxt,0))',Np,Np);
						time_Laplacian=(Nt*Dt*It*I_all)'*Mxt*(Nt*Dt*It*I_all);
						time_Laplacian=B2*inv_Mxt*B1T;
						for i=1:Nr
							temp(:)=0;
							temp(i)=1;
							v=P(temp);
							Sv=time_Laplacian*pr'*invA(Mxt*pr*v);
							S_tilde(:,i)=pr'*invA(Mxt*pr*Sv);
							S_tilde(:,i)=P(S_tilde(:,i));
						end

						inverse_S_tilde=sparse_inverse;
						ctrl_loc=ctrl_solver;
						index_agmg=1;
						ctrl_loc.init('direct',...
    											1e-13,...
    											ctrl_inner11.itermax,...
    											ctrl_inner11.omega,...
    											1,...
    											'S_tilde',...
													index_agmg);
						
						inverse_S_tilde.init(S_tilde+1e-10*speye(Nr,Nr),ctrl_loc);
						inverse_approx_S=@(y) inverse_S_tilde.apply(y);
					end
					
										%


					%
					schur_application = @(x) P(C*(x) + B2*(invA(B1T*x)));
					%schur_application = @(x) P(C*P(x) + B2*(invA(B1T*P(x))));
					%schur_application = @(x) P(C*P(x) + B2*(1.0./spdiags(A,0).*(B1T*P(x))));
					ctrl_loc=ctrl_solver;
					% passing index_agmg+1 we will use one copy of agmg, allowing preprocess 
					ctrl_loc.init('gmres',...
												1e-12,...
												1000,...
												0,...
												1,...
												'iter',...
												index_agmg);
					%inverse_approx_S = @(y) apply_iterative_solver( schur_application, y, ctrl_loc,@(z) z)% inverse_approx_S);

					invDiagA=sparse(1:Np,1:Np,(1./spdiags(A,0))',Np,Np);
					SCA = B2*invDiagA*B1T;
					inv_SCA=sparse_inverse;
					index_agmg=index_agmg+1;
					ctrl_loc=ctrl_solver;
					ctrl_loc.init('agmg',...
    										1e-6,...
    										100,...
    										0,...
    										1,...
    										'SCA',...
												index_agmg);
					inv_SCA.init(SCA + 0e-12*speye(Nr,Nr),ctrl_loc);
					inv_SCA.info_inverse.label='schur_ca';
					inv_SCA.cumulative_iter=0;
					inv_SCA.cumulative_cpu=0;
					

					%define inverse action
					%inverse_approx_S=	@(z) inv_SCA.apply(z);


					% gamma=0.00;
					% alpha=0;
					% comm=B2*inv_Mxt*B1T;
					% middle=gamma*C+B2*inv_Mxt*A*inv_Mxt*B1T;

					% inv_comm=sparse_inverse;
					% index_agmg=index_agmg+1;
					% ctrl_loc=ctrl_solver;
					% ctrl_loc.init('agmg',...
    			% 							1e-10,...
    			% 							100,...
    			% 							0,...
    			% 							1,...
    			% 							'SCA',...
					% 							index_agmg);
					% inv_SCA.init(comm + 0e-12*speye(Nr,Nr),ctrl_loc);
					% inv_SCA.info_inverse.label='schur_ca';
					% inv_SCA.cumulative_iter=0;
					% inv_SCA.cumulative_cpu=0;

					% inv_diag_SCA=sparse(1:Nr,1:Nr,(1./spdiags(SCA,0))',Nr,Nr);

					%inverse_approx_S=	@(z) inv_SCA.apply(middle*inv_SCA.apply(z))+alpha*diinv_diag_SCA*z;
					
					%assemb=prec_times_matrix(P,matrix_times_prec(SCA,P));
					%inverse_approx_S=	@(z) (assemb+1e-12*speye(Nr))\z;

					
					%
					% PREPROCESS NULL SPACE METHOD
					%
					time_null=tic;
					%inverse_approx_S=@(y) (S+1e-12*speye(Nr,Nr))\y;
					invS_Wmat=prec_times_matrix(inverse_approx_S,W_mat'); %S^{-1}Z
					MS=W_mat*(invS_Wmat);	 %MS=Z^TA^{-1} Z				
					inverse_MS=sparse_inverse; % create LU factorization of MS
					ctrl_loc=ctrl_solver;
					ctrl_loc.init('direct',...
					 							1e-12,...
					 							200,...
					 							0,...
					 							0,...
					 							'MS',...
					 							84);
					inverse_MS.init(MS,ctrl_loc);


					
					null_space_application=@(y) null_space_method(y,P,...
																												inverse_approx_S,...
																												W_mat,invS_Wmat,...
																												@(x) inverse_MS.apply(x));
					%inverse_S=@(y) -apply_Schur_inverse(y,null_space_application,Schur);
					inverse_S=@(y) -apply_Schur_inverse(y,inverse_approx_S,Schur);
					%inverse_S=@(y) -(Schur+1e-12)\y;


					%temp=matrix_times_prec(S,null_space_application);
					if (build_real_schur)
						%S=P(Stt+Sxx+Stx+Stx');
						temp=prec_times_matrix(inverse_S,-Schur);

						% %ragmgeturnagmg
						%figure
						[real,imag,eigenvec,perm]= study_eigenvalues(temp,'never gonna works');
						%plot(real,imag,'o')

						figure
						semilogy(abs(real))
						title('prec S')
						
						%plot(real,imag,'o')
						figure
						plot(real,imag,'o')

						real(Nr)/real(N+1)
						
					end


					%
					%	(C+Dt invA Dt^T)^{-1}=
					% (C+time_laplacian*invA)^{-1}=
					% (A_rho)*(C*A_rho+time_laplacian)^{-1}=:prec
					%

				

					% set the vectors in the kernels of A
					A_kernels=zeros(N+1,Np);
					ncellphi=ncell;
					for k=1:N+1
						A_kernels(k,1+(k-1)*ncellphi:k*ncellphi)=ones(1,ncellphi)/sqrt(ncellphi);
					end
					
					prec = @(x) SchurCA_based_preconditioner(x,@(y) ort_proj(invA(y),A_kernels),...
																									 @(z) inverse_S(z),...
																									 @(y) B1T*y,...
																									 @(z) P(B2*z),...
																									 Np,Nr,...
																									 'full',ncell);

					

					return
					J=[
						 A, B1T;
						 prec_times_matrix(P,B2), -prec_times_matrix(P,C)
					];

					rhsR=[f1;P(f2)];
				

					%sol=stationary_iterative_methods(@(x) J*x,rhsR,zeros(Np+Nr,1),1e-3,100,@(z) 0.1*prec(z));
					%res=J*sol-rhsR;
					%plot(res)
					%norm(res)/norm(rhsR)
					return

					
					temp=prec_times_matrix(prec,J);
					%temp=matrix_times_prec(J,prec);

				

					% %return
					%figure
					[real,imag,eigenvec,perm]= study_eigenvalues(temp,'never gonna works');
					%plot(real,imag,'o')

					figure
					semilogy(abs(real(1:N+1+N+2)))

					figure
					plot(real)

					
					
					figure
					plot(real,imag,'o')
					%semilogy(abs(real))
					title('no hope')
					fprintf('%1.1e, %1.1e (min,max) max/min %1.1e \n', real(N+1), real(Nr),	real(Nr)/		real(N+1))
				
											 
					return

					[real,imag,eigenvec,perm]= study_eigenvalues(my,'my');
					figure
					semilogy(abs(real(N+1:Nr)))
					title('my')
					fprintf('%1.1e, %1.1e (min,max) max/min %1.1e \n', real(N+1), real(Nr),	real(Nr)/		real(N+1))
					
					
					test=prec_times_matrix(@(y) inverse_SC.apply(y),my);

					[real,imag,eigenvec,perm]= study_eigenvalues(test,'prec_S');
					figure
					semilogy(abs(real(N+1:Nr)))
					title('stt/my')
					fprintf('%1.1e, %1.1e (min,max) max/min %1.1e \n', real(N+1), real(Nr),	real(Nr)/		real(N+1))
					return

					
					

					
					


					
					norm(full(Stt-my))

					figure
					plot(diag(Stt))
					figure
					plot(diag(my))

					
					return
					
					

					mode=1
					if (mode==1)
						invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);
						approx_SCA = sparse(C+B2*invDiagA*B1T);
						
						% passing index_agmg+1 we will use one copy of agmg, allowing preprocess
						inv_S=sparse_inverse;
						ctrl_loc=ctrl_solver;
						ctrl_loc.init('agmg',...
													1e-1,...
													200,...
													0,...
													0,...
													'C+B1diag(A)^{-1}B2',...
													73);
						inv_S.init(approx_SCA+1e-12*speye(Nr,Nr),ctrl_loc);
						inv_S.info_inverse.label='schur_ca';
						inverseS = @(x) inv_S.apply(x);
						%%
					elseif (mode==2)
						Pspace=assemble_space_projector(N+1,I');  
						Ptime=assemble_time_projector(N,ncell2h);
						P_pr = Pspace'*Ptime';
						
						time_Laplacian=(Nt*Dt*It*I_all)'*Mxt*(Nt*Dt*It*I_all);
						A_rho=Ptime*Pspace*A*Pspace'*Ptime';

						Stilde  = A_rho*C+time_Laplacian;

						inv_Stilde=sparse_inverse;
						ctrl_loc=ctrl_solver;
						ctrl_loc.init('agmg',...
													1e-4,...
													200,...
													0,...
													0,...
													'S_tilde',...
													73);
						inv_Stilde.init(Stilde+1e-12*speye(Nr,Nr),ctrl_loc);


						inverseS = @(x) inv_S.apply(x);


					end
						
					invS_Wmat=prec_times_matrix(inverseS,W_mat');
					MS=W_mat*(invS_Wmat);
					
					
					inverse_MS=sparse_inverse;
					ctrl_loc=ctrl_solver;
					ctrl_loc.init('direct',...
					 							1e-14,...
					 							200,...
					 							0,...
					 							0,...
					 							'MS',...
					 							84);
					inverse_MS.init(MS,ctrl_loc);

					prec_S=zeros(Nr,Nr);
					for i=1:Nr
						c=S_without_P(:,i);
						prec_S(:,i)=null_space_method(c,P,inverseS,W_mat,invS_Wmat,	inverse_MS);
					end

					[real,imag,eigenvec,perm]= study_eigenvalues(prec_S,'prec_S');
					figure
					semilogy(abs(real(N+1:Nr)))
					title('prec_null_space')
					fprintf('%1.1e, %1.1e (min,max) max/min %1.1e \n', real(N+1), real(Nr),	real(Nr)/		real(N+1))

					pure=prec_times_matrix(@(x) P(inverseS(x)),S);
					[real,imag,eigenvec,perm]= study_eigenvalues(pure,'prec_S');
					hold on
					semilogy(abs(real(N+1:Nr)))
					title('prec_S')
					
					[real,imag,eigenvec,perm]= study_eigenvalues(S,'prec_S');
					semilogy(abs(real(N+1:Nr)))
					title('S')
					hold off
					fprintf('%1.1e, %1.1e (min,max) max/min %1.1e \n', real(N+1), real(Nr),	real(Nr)/		real(N+1))
					
					fprintf(eigen_fid,'%d,%d, %1.1e, %1.1e \n',2^(h_i-1),2^(dt_i-1), real(N+1), real(Nr))
					return
					fclose(eigen_fid)
					
					
				end

				if (0)
					% time laplacian
					time_Laplacian=I_all'*Dt'*Nt*Nt*Dt*I_all;
					
					% first space, then time 
					Pspace=assemble_space_projector(N,I');  
					Ptime=assemble_time_projector(N,ncell);
					A_rho_ts=Pspace*Ptime*A*Ptime'*Pspace';
					
					% first time, then space
					Pspace=assemble_space_projector(N+1,I');  
					Ptime=assemble_time_projector(N,ncell2h);
					A_rho_st=Ptime*Pspace*A*Pspace'*Ptime';
					
					% A_rho_st and A_rho_ts are equal
 					A_rho=A_rho_st;
					% inverse of A_rho
					invA_rho=sparse_inverse;
					ctrl_loc=ctrl_solver;
					ctrl_loc.init('direct',...
												1e-14,...
												200,...
												0,...
												0,...
												'A_rho',100);
					invA_rho.init(A_rho+1e-12*speye(Nr,Nr),ctrl_loc);
								
				
					S_tilde=zeros(Nr,Nr);
					temp=zeros(Nr,1);
					for i=1:Nr
						temp(:)=0;
						temp(i)=1;
						v=P(temp);
						sum(v)
						v2=time_Laplacian*v;
						sum(v2)
						S_tilde(:,i)=P(invA_rho.apply(time_Laplacian*v)+C*v);
					end
				end

			

				
				
				

				if (1)

					% define the inverse 
					invA=sparse_inverse;
					ctrl_loc=ctrl_solver;
					ctrl_loc.init('agmg',...
												1e-14,...
												200,...
												0,...
												0,...
												'A',64);
					invA.init(A+1e-12*speye(Np,Np),ctrl_loc);
					
					
				

					mode=1
					if (mode ==1)
						time_Laplacian=(Nt*Dt*It*I_all)'*Mxt*(Nt*Dt*It*I_all);
					elseif (mode==2)
						time_Laplacian=(Nt*Dt*It*I_all)'*(Nt*Dt*It*I_all);
					elseif(mode==3)
						Dt_rho = assembleDt(N-2,ncell2h);
						time_Laplacian = (Nt*Dt_rho)'*(Nt*Dt_rho);
					end
					
					
					% define operator that project rho_variable into phi_variables
					Pspace=assemble_space_projector(N,I'); Ptime=assemble_time_projector(N,ncell);
					P_pr = (Pspace*Ptime)';
					
					%Pspace=assemble_space_projector(N+1,I'); Ptime=assemble_time_projector(N,ncell2h);
					%P_pr = Pspace'*Ptime';

					
					%
					S_tilde=zeros(Nr,Nr);
					temp=zeros(Nr,1);
					for i=1:Nr
						temp(:)=0;
						temp(i)=1;
						v=P(temp);
						if (mode ==1)
							S_tilde(:,i)=time_Laplacian*P_pr'*invA.apply(Mxt*P_pr*v);
							%S_tilde(:,i)=P_pr'*invA.apply(Mxt*P_pr*time_Laplacian*v);
						elseif(mode==2)
							%S_tilde(:,i)=time_Laplacian*P_pr'*Mxt*invA.apply(Mxt*P_pr*v);
							S_tilde(:,i)=P_pr'*Mxt*invA.apply(Mxt*P_pr*time_Laplacian*v);
						elseif(mode==3)
							S_tilde(:,i)=time_Laplacian*P_pr'*Mxt*invA.apply(Mxt*P_pr*v);
						end
							
						S_tilde(:,i)=P(S_tilde(:,i));
					end
				end

				%norm(Stt-S_tilde)/norm(Stt)

				
				% Schur complement S
				%S=PCP+Stt+Stx+Stx'+Sxx;
				S=Stt;
				
				% inverse of S
				invS=sparse_inverse;
				ctrl_loc=ctrl_solver;
				ctrl_loc.init('direct',...
											1e-14,...
											200,...
											0,...
											0,...
											'S',64);
				invS.init(S+1e-12*speye(Nr,Nr),ctrl_loc);
				inverseS=@(y) invS.apply(y);

				approx=S_tilde;
				
				my_S  = apply_prec_matrix(inverseS,approx);
				[real,imag,eigenvec]= study_eigenvalues(my_S,'my_S');
				fprintf('max/min %1.1e\n', real(Nr)/real(N+1));
				figure
				semilogy(real(N+1:Nr))
				title('my')
				return

				C_p_Stt=PCP+Stt;
				
				%my_S  = apply_prec_matrix(inverseS,C_p_Stt); see_eigen(my_S,'myS',N)
				
				
				C_S  = apply_prec_matrix(inverseS,PCP);
				Stt_S= apply_prec_matrix(inverseS,Stt);
				Stx_S= apply_prec_matrix(inverseS,Stx);
				Sxx_S= apply_prec_matrix(inverseS,Sxx);

				see_eigen(C_S,'CS',N)
				see_eigen(Stt_S,'Stt',N)
				see_eigen(Stx_S,'Stx',N)
				see_eigen(Sxx_S,'Sxx',N)
				
				


				
				return

				
				% assembly SAC=-(C+B2 * diag(A)^{-1} B1T)
				[prec] = assemble_prec_dual(A,B1T,B2,C,1e-1);
				
				precS=apply_prec_matrix(prec,S);
				
				[real,imag,eigenvec]= study_eigenvalues(precS,'prec_S');
				
				fprintf('max/min %1.1e\n', real(Nr)/real(N+1))
				figure
				plot(real,imag,'o')
				figure
				
			end 

			
			timelinsys=tic;
			[omegak, info_solver_newton,norm_ddd,resume_msg] = solvesys(JOC,OC, controls);
			%[res,resp,resr,ress]=compute_linear_system_residuum(JOC,OC,omegak);
			
			
			%state_message=sprintf('%d - | res p|=%1.4e |res r|=%1.4e |res s|=%1.4e |res|=%1.2e' ,...
			%			itk2+1, resp,resr,ress,res);
			%fprintf('%s \n',state_message);
			
			
			sum_linsys= sum_linsys+toc(timelinsys);
      sum_total = sum_total+toc(total);
			sum_prec  = sum_prec+info_solver_newton.prec_cpu;
			if (strcmp(resume_msg,''))
				print_info_solver(info_solver_newton)
				print_info_solver(info_solver_newton,logID)
			else
				fprintf('%s\n',resume_msg);
				fprintf(logID,'%s\n',resume_msg);
			end

			if (~ (info_solver_newton.flag == 0) )
				disp('ERROR')
				%return	
			end
			
			% deltat=1/(N+1);
			% masses=zeros(N,1);
			% for i = 1:N
			%   masses(i)=omegak(tnp+1+(i-1)*ncell2h:tnp+i*ncell2h)'*area2h;
			%   state_message=sprintf('y*area= %1.4e',min(masses(i)));
			%   fprintf('%s \n',state_message);
			%   fprintf(logID,'%s \n',state_message);
			% end
			
			%return
			

			
			sum_iter_outer_linear= sum_iter_outer_linear+info_solver_newton.outer_iter;
			sum_iter_inner_linear= sum_iter_inner_linear+info_solver_newton.inner_iter;
			sum_iter2_inner_linear= sum_iter2_inner_linear+info_solver_newton.inner_iter2;
			sum_iter3_inner_linear= sum_iter3_inner_linear+info_solver_newton.inner_iter3;


			
			
			% Linesearch just to ensure that rho and s stay positive
			alfak = 1;
			while any(uk(tnp+1:end)+alfak*omegak(tnp+1:end)<=0)
				alfak = 0.5*alfak;
				if alfak < alfamin
					alfak = alfamin;
					break
				end
			end

			uk = uk + alfak*omegak;
			
			%fprintf('Step=%1.2e\n',alfak)
			if alfak==alfamin
				disp('Step too small')
				flag2 = 1;
			end

			% se phi to sum to zero
			uk(1:Np)=ort_proj(uk(1:Np),ones(1,Np));
			
			% normalize rho
			for k = 1:N
				rhomass=area2h'*uk(Np+1+(k-1)*ncell2h:Np+k*ncell2h);
				%fprintf('mass rho(%d)-1.0=%1.1e\n',k,rhomass-1.0)
				uk(Np+1+(k-1)*ncell2h:Np+k*ncell2h) = uk(Np+1+(k-1)*ncell2h:Np+k*ncell2h)/rhomass;
			end

			if (save_data)
				data_name=sprintf('/DS%d_inner%d',itk1,itk2);
				h5create(filename_h5,data_name,[Np+2*Nr+2])
				h5write(filename_h5,data_name,[uk;mu;theta]')
			end

			
			%catch
			%  disp('Not solved, breaking')
			%  return
			%end
		end

    
    if flag_theta==1
      break
    end
    
    phimu = uk(1:tnp);
    rhomu = uk(tnp+1:tnp+tnr2h);
    smu = uk(tnp+tnr2h+1:end);

    
    
    % error bound on optimality
    delta_0 = (N/Nt)*sum(area)*mu;

    % 
    state_message=sprintf('%5s %3i %7s %1.4e %7s %3i  %10s %1.4e %4s %1.2e %7s %1.2e',...
	    'STEP:',itk1,...
	    ' Err',delta_0,...
	    ' NEWTON ',itk2,...
	    ' min(rho) ',min(uk(tnp+1:tnp+tnr2h)),' mu ',mu,' theta ',theta);
    if verb>0
        fprintf('%s \n', state_message)
	
    end
    cpu_total=cpu_total+sum_total;
    cpu_assembly=cpu_assembly+sum_assembly;
    cpu_linsys=cpu_linsys+sum_linsys;
    fprintf(logID,'%s \n',state_message);
  
    cost_message=sprintf('LINSYS: NEWTON %8.1e %2d OUT %3d IN %5d | CPU: LINSYS %1.4e ASSEMBLY %1.4e',...
			 delta_mu,itk2,sum_iter_outer_linear,uint64(sum_iter_inner_linear),...
			 sum_linsys,sum_assembly);
    fprintf('%s \n',cost_message);
    fprintf(logID,'%s \n',cost_message);

    fprintf(csvID,'%6d,%6d,%6d,%5d,%8.3e,%6d,%7d,%8d,%11d,%1.1e,%1.1e,%1.1e,%1.1e,%1.1e,%1.1e,%8d,%11d,%8d,%11d\n',...
	    ncell2h,ncell,N,...
	    itk1,delta_0,itk2,sum_iter_outer_linear,uint64(sum_iter_inner_linear),...
	   info_solver_newton.inner_nequ,min(uk(tnp+1:tnp+tnr2h)),mu,theta,...
	   sum_linsys,sum_assembly,sum_prec,...
	   uint64(sum_iter2_inner_linear),info_solver_newton.inner_nequ2,...
	   uint64(sum_iter3_inner_linear),info_solver_newton.inner_nequ3);

end

if (save_data)
  data_name=sprintf('/DS%d',itk1+1);
  h5create(filename_h5,data_name,[Np+2*Nr+2])
  h5write(filename_h5,data_name,[uk;mu;theta]')
end


fprintf('%17s %4i %7s %1.4e \n','Total Newton its:',tit,'Error: ',delta_0)
fprintf(controls.logID,'%17s %4i %7s %1.4e \n','Total Newton its:',tit,'Error: ',delta_0);


phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);

% Compute Wasserstein distance
rho_all = [rho_in;rho;rho_f];
W2th = compute_cost(ind,edges,mid,cc,gradt,Mst,RHt,It,N,rho_all,phi,rec);
fprintf(controls.logID,'%35s %1.4e \n','Approximated Wasserstein distance: ',W2th);
fprintf(controls.logID,'%19s %1.4e %21s %1.4e \n','Total linsys time: ',cpu_linsys,'Total assembly time: ',cpu_assembly);

if bc_sol==1&&compute_err==1
    
    % if the solution is known compute errors
    
    % shift the constant of phi
    [~,indc]=min((cc(:,1)-0.5).^2+(cc(:,2)-0.5).^2);
    phi = phi+(potential(0.5,0.5,0)-phi(indc));
    
    rhoa = RHt*It*rho_all;
    [rhos]=compute_rhosigma(ind,edges,cc,mid,N,rho_f,rho_in,gradt,Mst,RHt,It,Rst,rec,uk,'rhos');

    % Compute the errors
    err_cost = abs(W2-W2th); % cost error
    err_p = 0;
    err_rhoth = 0;
    for k=1:Nt
        t = (k-0.5)/Nt;
        phi_real = potential(cc(:,1),cc(:,2),t);
        rho_real = geodesic(cc2h(:,1),cc2h(:,2),t);
        rho_real = rho_real*mass/sum(rho_real.*area2h);
        rho_t = 0.5*(rho_all((k-1)*ncell2h+1:k*ncell2h)+rho_all(k*ncell2h+1:(k+1)*ncell2h));
        % error on the potential
        err_p = err_p + (1/Nt)*sum( rhoa((k-1)*ncell+1:k*ncell).*Mx*(phi_real-phi((k-1)*ncell+1:k*ncell)).^2 );
        % error on the geodesic
        err_rhoth = err_rhoth + (1/Nt)*sum(Mx2h*abs(rho_t-rho_real));
    end
    err_p = sqrt(err_p);
    
    fprintf(controls.logID,'%10s %1.4e \t %11s %1.4e \t %11s %1.4e \n','W2-error: ',err_cost,'phi-error: ',err_p,'rho-error: ',err_rhoth);
    
end

msg=sprintf('%10s %1.4e \t %11s %1.4e \t %11s %1.4e','W2-error: ',err_cost,'phi-error: ',err_p,'rho-error: ',err_rhoth);

fprintf(controls.logID,'%s\n',msg);
fprintf('%s\n',msg);
% plot
if (plot_figures)
figure
fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rho_in,'FaceColor','flat');
colorbar
axis square
axis off
caxis([0 max(rho)])
%str =num2str(0);
%outname=strcat('cross_harm_rhok',str,'.jpg');
%saveas(fig,outname)
figure
fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rho_f,'FaceColor','flat');
colorbar
axis square
axis off
caxis([0 max(rho)])
%str =num2str(nk+1);
%outname=strcat('cross_harm_rhok',str,'.jpg');
%saveas(fig,outname)
%dsp = 2;
dsp = ceil(N/2);
for k=dsp:dsp:N-dsp+1
    %k=4;
    rhok = rho((k-1)*ncell2h+1:k*ncell2h);
    figure
    fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rhok,'FaceColor','flat');
    colorbar
    axis square
    axis off
    caxis([0 max(rho)])
    %str =num2str(k);
    %outname=strcat('gauss_linear_rhok',str,'.jpg');
    %saveas(fig,outname)
end

end 


fclose(logID);
fclose(csvID);
%h5disp(filename_h5)
