%% Code starts
if mesh_type == 1
    mesh_name = strcat('meshes/tri2_mesh',num2str(h_i));
elseif mesh_type == 2
    mesh_name = strcat('meshes/subtri2_mesh',num2str(h_i));
else
    mesh_name = strcat('meshes/sq_mesh',num2str(h_i));
end

cpu_total=0.0;
cpu_assembly=0.0;


% load mesh
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

filename=strcat('runs/',str_test,approach_string);%,controls_string);
disp(filename)
log_filename=strcat(filename,'.log')
logID = fopen(log_filename,'w');
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
fprintf(csvID,'  nrho,    np,    nt, step,    error,newton,  outer,   inner,  innernequ,  minrho,      mu,   theta, cpulinsys, cpuassemb\n');

% Boundary conditions:
[rho_in,rho_f,mass] = bc_density(test_case,cc2h,area2h);

% Barrier method's parameters:
eps_0 = 1e-6; % tolerance for the final solution
k2max = 30; % maximum number of inner (Newton) iterations
k1max = 20; % maximum number of outer iterations
theta0 = 0.2; % decay ratio for the perturbation parameter mu
alfamin = 0.1; % minimal step length accepted


tnp = (N+1)*ncell; % total number of dofs for the potential
tnr2h = N*ncell2h; % total number of dofs for the density

mu0 = 1;
phi0 = zeros(tnp,1);
rho0 = (mass/sum(area))*ones(tnr2h,1);
s0 = mu0./rho0;
uk = [phi0;rho0;s0]; % initial condition of the barrier method

save_data=controls.save_data;
if (save_data)
  filename_save=strcat('runs/','PhiRhoSMuTheta',str_test,approach_string,'.dat');%,controls_string);
  IDsave = fopen(filename_save,'w');
  IDsave = write2td_sequence( IDsave, [uk;mu0;theta0], 0.0,'head');

  %filename_save=([strcat('runs/','PhiRhoSMuTheta',str_test,approach_string,'.dat')]);%,controls_string);
  %matvar = matfile(filename_save,'Writable',true);
  %save(filename_save,'matvar.[uk;mu0;theta0]','-append'); 
end
% Assemble matrices

Mx = spdiags(area,0,ncell,ncell);
Mx2h = spdiags(area2h,0,ncell2h,ncell2h);
ds = edges(ind.internal,5).*edges(ind.internal,6);
Ms = spdiags(ds,0,nei,nei);
div = Div2D(ncell,nei,ind,edges); % divergence matrix
grad = -Ms\div'; % gradient matrix


% global matrices
Dt = assembleDt(N,ncell);
divt = assembleDivt(N,ncell,nei,div);
gradt = assembleGradt(N,ncell,nei,grad);
Mxt = assembleMxt(N,ncell,Mx);
Mxt2h = assembleMxt(N,ncell2h,Mx2h);
Mst = assembleMst(N,nei,Ms);
RHt = assembleRHt(N,ncell);
It = assembleIt(N,ncell,ncell2h,I);


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
        mu0 = 5*mu0;
        theta = theta0;
        s0 = mu0./rho0;
        mu = mu0/theta;
        itk1 = 0;
        tit = 0;
        uk = [phi0;rho0;s0];
        delta_0 = 2*eps_0;
        continue
    end
    
    mu = theta*mu;
    eps_mu = eps_0;
    itk2 = 0;
    flag2 = 0;

    sum_assembly=0;
    sum_total=0;
    sum_linsys=0;
    sum_iter_newton=0;
    sum_iter_outer_linear=0;
    sum_iter_inner_linear=0;
    while true
        total=tic;
        assembly=tic;
	
        
        OC = Fkgeod(ind,edges,cc,mid,N,(rho_f+mu)/(1+mu),(rho_in+mu)/(1+mu),Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk,mu);
        delta_mu = norm([OC.p;OC.r;OC.s]);

	state_message=sprintf('%d - |OC.p|=%1.4e |OC.r|=%1.4e |OC.s|=%1.4e \n' ,...
				itk2+1, norm(OC.p),norm(OC.r),norm(OC.s));
	fprintf(logID,'%s \n',state_message);
	if verb>1
	  if (  itk2 == 0)
	    state_message=sprintf(' \n')	;		
	    fprintf('%s \n',state_message);
	  end	
	  fprintf('%s \n',state_message);
	end



	% for i = 1:N
	%   masses(i)=uk(tnp+1+(i-1)*ncell2h:tnp+i*ncell2h)'*area2h;
	% end
	% state_message=sprintf('%1.4e<=masses rho[:]<= %1.4e',min(masses),max( masses));
	% fprintf('%s \n',state_message);
	% fprintf(logID,'%s \n',state_message);

	
	

	% current_res_phi = OC.p;
	% for i = 1:Nt
	%   imbalance_res_phi(i)=sum(current_res_phi((i-1)*ncell+1:i*ncell));%/norm(current_res_phi((i-1)*ncell+1:i*ncell));
	% end
	% state_message=sprintf('%1.4e<=sum(F_phi[:])/norm(F_phi[:])<= %1.4e',min(imbalance_res_phi),max(imbalance_res_phi));
	% fprintf('%s \n',state_message);
	% fprintf(logID,'%s \n',state_message);


	% current_res_phi = OC.r;
	% for i = 1:N
	%   imbalance_res_phi(i)=current_res_phi((i-1)*ncell2h+1:i*ncell2h)'*area2h;
	% end
	% state_message=sprintf('%1.4e<=integral_F_rho[:]<= %1.4e',min(imbalance_res_phi),max(imbalance_res_phi));
	% fprintf('%s \n',state_message);
	% fprintf(logID,'%s \n',state_message);

	

        if verb>1
            fprintf('%11s %3i %7s %1.4e \n','Inner step:',itk2,' Error:',delta_mu)
        end

        if delta_mu < eps_mu  
            if itk2<4
                theta = 0.8*theta;
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
            mu = mu/theta;
            theta = theta+0.2*(1-theta);
            mu = theta*mu;
            uk = [phimu;rhomu;smu];
            tit = tit+itk2;
            itk2 = 0;
            flag2 = 0;
            continue
        end

        % Compute the jacobian of the system of equations
        JOC = JFkgeod(ind,edges,cc,mid,N,...
		      (rho_f+mu)/(1+mu),(rho_in+mu)/(1+mu),...
		      Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk);
	sum_assembly=sum_assembly+toc(assembly);

	
	% Solve the linear system
	timelinsys=tic;
        [omegak, info_solver_newton] = solvesys(JOC,OC, controls,logID);
	sum_linsys= sum_linsys+toc(timelinsys);
        sum_total = sum_total+toc(total);
	print_info_solver(info_solver_newton)
	print_info_solver(info_solver_newton,logID)

	if (restart)
	  solver_approach=11;

	  ctrl_inner11=ctrl_solver;
	  ctrl_inner22=ctrl_solver;
	  ctrl_outer=ctrl_solver;
	  
          % krylov based solvers for M [x;y] = [f;g]
 	  % pcg works only with full
	  outer_solvers={'bicgstab'  ,'gmres','fgmres' ,'pcg'};

      % set here fgmres (for non stationary prec), bicgstab,gmres, pcg
	  for isolver=[3]%1:length(outer_solvers)
	    ctrl_outer.init(outer_solvers{isolver},controls.ctrl_outer.tolerance,3000,0.0,0); % verbose=[1,2] works only for fgmres
	    %ctrl_outer.init(outer_solvers{isolver},1e-10,3000,0.0,0); % verbose=[1,2] works only for fgmres
	    
	    
				% external prec appraoch
	    outer_precs={'full' ,'lower_triang'  ,'upper_triang','identity'};
	    nouter_precs=length(outer_precs);
	    
				% pcg works only with full
	    if (strcmp(outer_solvers,'pcg'))
	      nouter_precs=1
	    end
	    
				% fgmres works only with right prec
				%if (strcmp(outer_solvers,'fgmres'))
				%  nouter_precs=1
				%end
	    left_right='right';
	    
	    
	    for iprec=[3]%1:nouter_precs
	      outer_prec=outer_precs{iprec};
	      
		       % set here other approximate inverse of block11
	      ctrl_inner11.init('agmg',... %approach
				1e-12,... %tolerance
				10,...% itermax
				0.0,... %omega
				0,... %verbose
				'agmg10'); %label
					 %extra_info='full';
	      extra_info='block';
	      relax4_inv11=1e-12;
	      
	% set grounded_node>0 to gorund the potential in grounded node
	      indc=0;
	      indeces=set_grounding_node(JF.pp,ncellphi);

	      % ground other potential
	      grounded_values=omegak(indeces);
	      		      % set here list of solvers for block 22 
	      solvers={'agmg' ,'agmg'  ,'agmg' ,'direct','krylov' ,'krylov'  ,'incomplete','diag'};
	      iters  ={1      ,10      ,100,1       ,1        ,10        ,  1          ,0};
	      label  ={'agmg1','agmg10','agmg100','direct','krylov1','krylov10','incomplete','diag'};
	      relax4_inv22=0;
	      
	      for i=[2];%1:length(solvers)
		ctrl_inner22.init(solvers{i},1e-13,iters{i},1.0,0,label{i});
		restart_controls = struct('indc',grounded_node,...
				  'grounded_values',grounded_values,...
				  'sol',solver_approach,...
				  'outer_prec',outer_prec,...
				  'left_right',left_right,...
				  'ctrl_inner11',ctrl_inner11,...
				  'ctrl_inner22',ctrl_inner22,...
				  'ctrl_outer',ctrl_outer,...
				  'compute_eigen',compute_eigen,...
				  'verbose',verbose,...
				  'extra_info',extra_info,...
				  'relax4inv11',relax4_inv11,...
				  'relax4inv22',relax4_inv22);


	      end
	    end
	  end
	  timelinsys=tic;
	  omegak_old=omegak;
          [omegak, info_solver_newton] = solvesys(JOC,OC, restart_controls,logID);
	  sum_linsys=sum_linsys+toc(timelinsys);
	   print_info_solver(info_solver_newton);
	   print_info_solver(info_solver_newton,logID);
	  
	  if( 0 )
	  for i=1:Nt
	    grounded_values(i)=omegak(grounded_node+(i-1)*ncell);
	   fprintf('diff dir %1.4e diff phi%1.4e\n', ...
		    abs(restart_controls.grounded_values(i)-omegak(grounded_node+(i-1)*ncell)),...
		    norm(omegak(1+(i-1)*ncell:i*ncell)-omegak_old(1+(i-1)*ncell:i*ncell))...
		   )
	  end
	  for i=1:N
	    grounded_values(i)=omegak(grounded_node+(i-1)*ncell);
	    fprintf('diff rho %1.4e\n', ...
		    norm(omegak(Np+1+(i-1)*ncell2h:Np+i*ncell2h)-omegak_old(Np+1+(i-1)*ncell2h:Np+i*ncell2h))...
		   )
	  end
	  end
	  
	  %i=1;
	  %X=linspace(0,1,ncell)';
	  %diff=-omegak(controls.indc)+omegak(1+(i-1)*ncell:i*ncell)-omegak_old(1+(i-1)*ncell:i*ncell);
	  %plot(X,diff)

	  
	  fprintf('diff  sol = %1.4e \n', norm(omegak-omegak_old));
	  fprintf(logID,'diff  sol = %1.4e \n', norm(omegak-omegak_old));

	  
	 
	  
	end



	

	
	sum_iter_outer_linear= sum_iter_outer_linear+info_solver_newton.outer_iter;
	sum_iter_inner_linear= sum_iter_inner_linear+info_solver_newton.inner_iter;
	


	% y=omegak(Np+1:Np+Nr);
	% B1T = sparse(JOC.pr);
	% v=B1T*y;
	% for i = 1:Nt
	%   masses_increment(i)=sum(v(1+(i-1)*ncell:i*ncell));
	% end

	% state_message=sprintf('%1.4e<=sum(B1T y_i)[:]<= %1.4e',min(masses_increment),max( masses_increment));
	% fprintf('%s \n',state_message);
	% fprintf(logID,'%s \n',state_message);

	% v=OC.p-v;
	% for i = 1:Nt
	%   masses_increment(i)=sum(v(1+(i-1)*ncell:i*ncell));
	% end
	% state_message=sprintf('%1.4e<=sum(f-B1T y_i)[:]<= %1.4e',min(masses_increment),max( masses_increment));
	% fprintf('%s \n',state_message);
	% fprintf(logID,'%s \n',state_message);

	

	% for i = 1:N
	%   masses_increment(i)=omegak(tnp+1+(i-1)*ncell2h:tnp+i*ncell2h)'*area2h;
	% end

	% state_message=sprintf('%1.4e<=masses increment rho[:]<= %1.4e',min(masses_increment),max( masses_increment));
	% fprintf('%s \n',state_message);
	% fprintf(logID,'%s \n',state_message);


	%OC = Fkgeod(ind,edges,cc,mid,N,(rho_f+mu)/(1+mu),(rho_in+mu)/(1+mu),Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk,mu);
        %delta_mu = norm([OC.p;OC.r;OC.s]);


	
	
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
       
        if alfak==alfamin
            flag2 = 1;
        end

    end

    phimu = uk(1:tnp);
    rhomu = uk(tnp+1:tnp+tnr2h);
    smu = uk(tnp+tnr2h+1:end);

    if (save_data)
      IDsave=write2td_sequence( IDsave, [uk;mu;theta], itk1,'body');
      %s = size(m, 'my_variables');
      %matvar.my_variables(s+1, :) = [uk;mu;theta];
    end
    
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
    fprintf(logID,'%s \n',state_message);
  
    cost_message=sprintf('LINSYS: NEWTON %8.1e %2d OUT %3d IN %5d | CPU: LINSYS %1.4e ASSEMBLY %1.4e',...
			 delta_mu,itk2,sum_iter_outer_linear,uint64(sum_iter_inner_linear),...
			 sum_linsys,sum_assembly);
    fprintf('%s \n',cost_message);
    fprintf(logID,'%s \n',cost_message);

    fprintf(csvID,'%6d,%6d,%6d,%5d,%8.3e,%6d,%7d,%8d,%11d,%1.2e,%1.2e,%1.2e,%1.4e,%1.4e \n',...
	    ncell2h,ncell,N,...
	    itk1,delta_0,itk2,sum_iter_outer_linear,uint64(sum_iter_inner_linear),...
	   info_solver_newton.inner_nequ,min(uk(tnp+1:tnp+tnr2h)),mu,theta,...
	   sum_linsys,sum_assembly);

end

if (save_data)
  IDsave=write2td_sequence( IDsave, [uk;mu;theta], itk2,'tail');
%  fclose(IDsave);
end


fprintf('%17s %4i %7s %1.4e \n','Total Newton its:',tit,'Error: ',delta_0)

phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);

% Compute Wasserstein distance
rho_all = [rho_in;rho;rho_f];
W2 = compute_cost(ind,edges,mid,cc,gradt,Mst,RHt,It,N,rho_all,phi,rec);

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


