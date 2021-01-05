clear
close all


% Convergence test
% test 1 -> pure translation
% test 2 -> contraction
test = 2;

% Number of refinement steps
nstep = 2;
% Initial number of intermediate densities
N_0 = 1;
% Initial level of mesh refinement
i_h = 1;

% Type of mesh
mesh_type = 1;

% Type of reconstruction
rec = 1;

% Verbosity level: {0,1,2}
verb = 0;


%% Code starts

if test==1
    supp = pi/0.3;
    initial =@(x,y) (((x-0.3).^2+(y-0.3).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.3).^2+(y-0.3).^2))));
    final =@(x,y) (((x-0.7).^2+(y-0.7).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.7).^2+(y-0.7).^2))));
    midpoint = @(x,y) (((x-0.5).^2+(y-0.5).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.5).^2+(y-0.5).^2))));
    potential=@(x,y,t) 0.4*x+0.4*y-(0.4^2)*t-0.4;
    geodesic =@(x,y,t) (((x-0.3-t*0.4).^2+(y-0.3-t*0.4).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.3-t*0.4).^2+(y-0.3-t*0.4).^2))));
    W2cost = sqrt(2)*0.4;
else
    sc = 0.3;
    initial =@(x,y) (1+cos(2*pi*(x-0.5)));
    final =@(x,y) ((((x-0.5)/sc).^2-0.5^2)<=0).*(1+cos(2*pi*((x-0.5)/sc)))/sc;
    midpoint =@(x,y) ((((x-0.5)/(0.5*sc+0.5)).^2-0.5^2)<=0).*(1+cos(2*pi*((x-0.5)/(0.5*sc+0.5))))/(0.5*sc+0.5);
    potential =@(x,y,t)((sc-1)/((sc-1)*t+1))*0.5*(x-0.5).^2;
    geodesic =@(x,y,t) ((abs(x-0.5)-0.5*(t*(sc-1)+1))<=0).*(1+cos(2*pi*((x-0.5)/(t*(sc-1)+1))))/(t*(sc-1)+1);
    W2cost = sqrt(((pi^2-6)*(sc-1)^2)/(12*pi^2));
end

if mesh_type == 1
    mesh_name = 'meshes/tri2_mesh';
elseif mesh_type == 2
    mesh_name = 'meshes/subtri2_mesh';
else
    mesh_name = 'meshes/sq_mesh';
end

% Barrier method's parameters:
eps_0 = 1e-8; % tolerance for the final solution
k2max = 30; % maximum number of inner (Newton) iterations
k1max = 20; % maximum number of outer iterations
theta0 = 0.2; % decay ratio for the perturbation parameter mu
alfamin = 0.1; % minimal step length accepted

hstep = zeros(nstep,1);
hstep2h = zeros(nstep,1);
Nstep = zeros(nstep,1);
cost = zeros(nstep,1);
err_cost = zeros(nstep,1);
err_mid = zeros(nstep,1);
err_rhoth = zeros(nstep,1);
err_p = zeros(nstep,1);
err_gp = zeros(nstep,1);

for i=1:nstep

    % load mesh
    fn = strcat(mesh_name,num2str(i_h+i-1));
    load(fn)
    
    if exist('S','var')==0
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

    
    % find the closest node to (0.5,0.5) to use to impose the condition on
    % the potential
    [~,indc]=min((cc(:,1)-0.5).^2+(cc(:,2)-0.5).^2);
    
    
    rho_in = initial(cc2h(:,1),cc2h(:,2));
    mass = sum(area2h.*rho_in);
    rho_in = rho_in/mass; mass = sum(area2h.*rho_in);
    rho_f = final(cc2h(:,1),cc2h(:,2));
    rho_f = rho_f*mass/sum(rho_f.*area2h);
    rho_m = midpoint(cc2h(:,1),cc2h(:,2));
    rho_m = rho_m*mass/sum(rho_m.*area2h);
    
    N = (N_0+1)*2^(i-1)-1; %number of intermediate points, >=1, odd
    Nt = N+1;
    
    hstep(i) = h;
    hstep2h(i) = h2h;
    Nstep(i) = N+1;
    
    tnp = (N+1)*ncell; % total number of dofs for the potential
    tnr2h = N*ncell2h; % total number of dofs for the density

    mu0 = 1;
    phi0 = zeros((N+1)*ncell,1);
    rho0 = (mass/sum(area))*ones(N*ncell2h,1);
    s0 = mu0./rho0;
    uk = [phi0;rho0;s0];

    
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

    while true

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

        while true

            OC = Fkgeod(ind,edges,cc,mid,N,(rho_f+mu)/(1+mu),(rho_in+mu)/(1+mu),Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk,mu);
            delta_mu = norm([OC.p;OC.r;OC.s]);

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

            % Solve the linear system
            omegak = solvesys(JOC,OC,indc);

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

        if verb>0
            fprintf('%11s %3i %7s %1.4e %9s %3i %10s %1.4e %4s %1.2e %7s %1.2e \n','OUTER STEP:',itk1,' Error:',...
                delta_0,' #NewIts:',itk2,' min(rho):',min(uk(tnp+1:tnp+tnr2h)),' mu:',mu,' theta:',theta)
        end

    end

    fprintf('%17s %4i %7s %1.4e \n','Total Newton its:',tit,'Error: ',delta_0)

    phi = uk(1:tnp);
    rho = uk(tnp+1:tnp+tnr2h);
    rhomid = rho((N-1)/2*ncell2h+1:(N-1)/2*ncell2h+ncell2h);

    % Compute the Wasserstein distance
    rho_all = [rho_in;rho;rho_f];
    gradphi = gradt*phi;
    if rec==1
        Rst = sparse((N+1)*nei,(N+1)*ncell);
        for k=1:N+1
            Rst((k-1)*nei+1:k*nei,(k-1)*ncell+1:k*ncell) = Ktos2D(ind,edges,cc,mid);
        end
        rhoa = RHt*It*rho_all;
        rhos = Rst*rhoa;
    else
        rhoa = RHt*It*rho_all;
        rhos = zeros((N+1)*nei,1);
        for k=1:N+1
            rhos((k-1)*nei+1:k*nei) = rho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
        end
    end
    cost(i) = (1/Nt)*rhos'*Mst*(gradphi.^2); cost(i) = sqrt(cost(i));

    % Compute the errors
    err_cost(i) = abs(W2cost-cost(i)); % cost error
    err_mid(i) = sum(Mx2h*abs(rhomid-rho_m)); % density error in the midpoint
    for k=1:Nt
        t = (k-0.5)/Nt;
        phi_real = potential(cc(:,1),cc(:,2),t);
        gphi_real = grad*phi_real;
        rho_real = geodesic(cc2h(:,1),cc2h(:,2),t);
        rho_real = rho_real*mass/sum(rho_real.*area2h);
        rho_t = 0.5*(rho_all((k-1)*ncell2h+1:k*ncell2h)+rho_all(k*ncell2h+1:(k+1)*ncell2h));
        % error on the potential
        err_p(i) = err_p(i) + (1/Nt)*sum( rhoa((k-1)*ncell+1:k*ncell).*Mx*(phi_real-phi((k-1)*ncell+1:k*ncell)).^2 );
        % error on the gradient of the potential
        err_gp(i) = err_gp(i) + (1/Nt)*sum( rhos((k-1)*nei+1:k*nei).*Ms*(gphi_real-gradphi((k-1)*nei+1:k*nei)).^2);
        % error on the geodesic
        err_rhoth(i) = err_rhoth(i) + (1/Nt)*sum(Mx2h*abs(rho_t-rho_real));
    end
    err_p(i) = sqrt(err_p(i));
    err_gp(i) = sqrt(err_gp(i));
    
    
    % print mid point
    fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rhomid,'FaceColor','flat');
    colorbar
    axis equal
    axis off
    str =num2str(i);
    outname=strcat('mid',str,'.jpg');
    saveas(fig,outname)
    close all
    

    clear I
    
end

A = [hstep2h hstep Nstep cost err_cost err_p err_gp err_mid err_rhoth];
rates = log(A(1:end-1,5)./A(2:end,5))./log(A(1:end-1,1)./A(2:end,1));
rates = [rates log(A(2:end,6)./A(1:end-1,6))./log(A(2:end,1)./A(1:end-1,1))];
rates = [rates log(A(2:end,7)./A(1:end-1,7))./log(A(2:end,1)./A(1:end-1,1))];
rates = [rates log(A(2:end,8)./A(1:end-1,8))./log(A(2:end,1)./A(1:end-1,1))];
rates = [rates log(A(2:end,9)./A(1:end-1,9))./log(A(2:end,1)./A(1:end-1,1))];
rates = [[0 0 0 0 0]; rates];
A = [A(:,1:5) rates(:,1) A(:,6) rates(:,2) A(:,7) rates(:,3) A(:,8) rates(:,4) A(:,9) rates(:,5)];

%save('A_tri2_comp_linear','A')


