clear
close all

% Computation of Optimal Transport with finite volumes
% (https://hal.archives-ouvertes.fr/hal-03032446)

% Time discretization. Uniform grid with Nt=N+1 intervals. N is the number
% of intermediate densities. Choose N odd, N>=1, in order to have the
% approximate geodesic midpoint.
N = 1;
Nt = N+1;

% Space discretization. Three types of mesh families available:
% 1 -> regular triangulation of the domain, with only acute angles
%      (https://www.i2m.univ-amu.fr/fvca5/benchmark/Meshes/index.html)
% 2 -> nested-constructed meshes, based on the previous family
%      (find details at https://hal.archives-ouvertes.fr/hal-03032446)
% 3 -> cartesian grids
% For each mesh, five levels of refinement h_i, 1->5, are available.
mesh_type = 1;
h_i = 1;
% Mesh structure:
% nodes -> array of nodes coordinates [x y]
% cells -> array of cells nodes [#nodes node1 node2 node3 ...]
% edges -> array of edges [node1 node2 K L d_sigma m_sigma m_sigma/d_sigma]
% ind -> indices struct: ind.internal -> indices of internal edges
%                        ind.bound    -> indices of boundary edges 
% cc -> array of cell centers coordinates [x y]
% area -> array of cell measures
% mid -> array of midpoints coordinates [x y]
% h -> meshsize, max(diam(K)) or max(sqrt(area))
% The outer mesh structure, if any, is labelled with 2h to distinguish it
% from the inner nested one


% Type of reconstruction for the density on the diamond cells.
% 1 -> weighted arithmetic mean
% 2 -> weighted harmonic mean
rec = 1;

verb = 1; % verbosity level: {0,1,2}


%% Code starts

if mesh_type == 1
    mesh_name = strcat('../meshes/tri2_mesh',num2str(h_i));
elseif mesh_type == 2
    mesh_name = strcat('../meshes/subtri2_mesh',num2str(h_i));
else
    mesh_name = strcat('../meshes/sq_mesh',num2str(h_i));
end

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

% Boundary conditions:
[rho_in,rho_f,mass] = bc_density(cc2h,area2h);

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
        omegak = solvesys(JOC,OC,1);
        
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

% Compute Wasserstein distance
rho_all = [rho_in;rho;rho_f];
W2 = compute_cost(ind,edges,mid,cc,gradt,Mst,RHt,It,N,rho_all,phi,rec);

% plot 
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





