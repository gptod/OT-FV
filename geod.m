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
fprintf(logID,'solver approach = %d\n',solver_approach);
ctrl_inner11.info(logID);
ctrl_inner22.info(logID);
ctrl_outer.info(logID);


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



	for i = 1:N
	  masses(i)=uk(tnp+1+(i-1)*ncell2h:tnp+i*ncell2h)'*area2h;
	end
	state_message=sprintf('%1.4e<=masses rho[:]<= %1.4e',min(masses),max( masses));
	%fprintf('%s \n',state_message);
	fprintf(logID,'%s \n',state_message);
	

	current_res_phi = OC.p;
	for i = 1:Nt
	  imbalance_res_phi(i)=sum(current_res_phi((i-1)*ncell+1:i*ncell));
	end
	state_message=sprintf('%1.4e<=imbalance_F_phi[:]<= %1.4e',min(imbalance_res_phi),max(imbalance_res_phi));
	fprintf('%s \n',state_message);
	fprintf(logID,'%s \n',state_message);

	for i = 1:Nt
	  imbalance_res_phi(i)=norm(current_res_phi((i-1)*ncell+1:i*ncell));
	  %state_message=sprintf('norm_F_phi( %d )= %1.4e',i,imbalance_res_phi(i));
	  %fprintf('%s \n',state_message);
	  %fprintf(logID,'%s \n',state_message);
	end
	state_message=sprintf('%1.4e<=norm_F_phi[:]<= %1.4e',min(imbalance_res_phi),max(imbalance_res_phi));
	fprintf('%s \n',state_message);
	fprintf(logID,'%s \n',state_message);

	

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
        JOC = JFkgeod(ind,edges,cc,mid,N,(rho_f+mu)/(1+mu),(rho_in+mu)/(1+mu),Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk);
	sum_assembly=sum_assembly+toc(assembly);

	
				% Solve the linear system
	timelinsys=tic;


	
	

	
        [omegak, info_solver] = solvesys(JOC,OC,controls,logID);
	sum_linsys=sum_linsys+toc(timelinsys);
	sum_total=sum_total+toc(total);
	sum_iter_outer_linear= sum_iter_outer_linear+info_solver.outer_iter;
	sum_iter_inner_linear= sum_iter_inner_linear+info_solver.inner_iter;
	
	print_info_solver(info_solver)
	print_info_solver(info_solver,logID)

	for i = 1:N
	  masses_increment(i)=omegak(tnp+1+(i-1)*ncell2h:tnp+i*ncell2h)'*area2h;
	end

	state_message=sprintf('%1.4e<=masses increment rho[:]<= %1.4e',min(masses_increment),max( masses_increment));
	%fprintf('%s \n',state_message);
	fprintf(logID,'%s \n',state_message);


	OC = Fkgeod(ind,edges,cc,mid,N,(rho_f+mu)/(1+mu),(rho_in+mu)/(1+mu),Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk,mu);
        delta_mu = norm([OC.p;OC.r;OC.s]);


	
	
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
	   info_solver.inner_nequ,min(uk(tnp+1:tnp+tnr2h)),mu,theta,...
	   sum_linsys,sum_assembly);

end

fprintf('%17s %4i %7s %1.4e \n','Total Newton its:',tit,'Error: ',delta_0)

phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);

% Compute Wasserstein distance
rho_all = [rho_in;rho;rho_f];
W2 = compute_cost(ind,edges,mid,cc,gradt,Mst,RHt,It,N,rho_all,phi,rec);

% plot
if (plot)
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


