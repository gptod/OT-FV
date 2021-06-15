clear
close all


% Convergence test
% test 1: sin -> pure translation
% test 2: compression -> contraction
test_case='sin';

% Type of mesh
mesh_type = 1;

% Type of reconstruction
% 1 -> weighted arithmetic mean
% 2 -> weighted harmonic mean
rec = 1;

% number of refinement
nstep = 1;

eps_0 = 1e-8; % tolerance 

verb = 0;

plot_figures=0
restart=0;
save_data=0;
compute_eigen=0;

verbose=0;




if strcmp(test_case,'sin')
    supp = pi/0.3;
    midpoint = @(x,y) (((x-0.5).^2+(y-0.5).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.5).^2+(y-0.5).^2))));
    potential=@(x,y,t) 0.4*x+0.4*y-(0.4^2)*t-0.4;
    geodesic =@(x,y,t) (((x-0.3-t*0.4).^2+(y-0.3-t*0.4).^2-0.3^2)<=0).*(1+cos(supp*(sqrt((x-0.3-t*0.4).^2+(y-0.3-t*0.4).^2))));
    W2 = sqrt(2)*0.4;
elseif strcmp(test_case,'compression')
    sc = 0.3;
    midpoint =@(x,y) ((((x-0.5)/(0.5*sc+0.5)).^2-0.5^2)<=0).*(1+cos(2*pi*((x-0.5)/(0.5*sc+0.5))))/(0.5*sc+0.5);
    potential =@(x,y,t)((sc-1)/((sc-1)*t+1))*0.5*(x-0.5).^2;
    geodesic =@(x,y,t) ((abs(x-0.5)-0.5*(t*(sc-1)+1))<=0).*(1+cos(2*pi*((x-0.5)/(t*(sc-1)+1))))/(t*(sc-1)+1);
    W2 = sqrt(((pi^2-6)*(sc-1)^2)/(12*pi^2));
else
    disp('error: wrong test case')
end


hstep = zeros(nstep,1);
hstep2h = zeros(nstep,1);
Nstep = zeros(nstep,1);
cost = zeros(nstep,1);
err_cost = zeros(nstep,1);
err_mid = zeros(nstep,1);
err_rhoth = zeros(nstep,1);
err_p = zeros(nstep,1);
err_gp = zeros(nstep,1);

for h_i = 1:nstep
    
    N = 2*2^(h_i-1)-1;
    %N = 4*2^(h_i-1)-1;
    Nt = N+1;

    str_test=sprintf('%s_h%d_rec%d_N%0.5d_',test_case,h_i,rec,N)


    solver_approach=11;

	ctrl_inner11=ctrl_solver;
	ctrl_inner22=ctrl_solver;
	ctrl_outer=ctrl_solver;
    
    ctrl_outer.init('fgmres',1e-5,3000,0.0,0); % verbose=[1,2] works only for fgmres
    left_right='right';

    outer_prec='full';

    % set here other approximate inverse of block11
    ctrl_inner11.init('diag',... %approach
            1e-12,... %tolerance
            10,...% itermax
            0.0,... %omega
            0,... %verbose
            'diag'); %label
    %extra_info='full';
    extra_info='block';
    relax4_inv11=1e-12;

    % set grounded_node>0 to gorund the potential in grounded node
    grounded_node=0;

    relax4_inv22=0;

    ctrl_inner22.init('agmg',1e-13,1,1.0,0,'agmg1');
    controls = struct('save_data',save_data,...
              'indc',grounded_node,...
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

    approach_string=strcat('schurCAwithdiagA_',...
               ctrl_outer.approach,'_',...
               left_right,'_',outer_prec,'_prec_',...
               'invA',ctrl_inner11.label,'_',...
               'invSCA',ctrl_inner22.label);

    disp(approach_string)



    geod;
    
    hstep(h_i) = h;
    hstep2h(h_i) = h2h;
    Nstep(h_i) = N;
    
    gradphi = gradt*phi;
    rhomid = rho((N-1)/2*ncell2h+1:(N-1)/2*ncell2h+ncell2h);
    rho_m = midpoint(cc2h(:,1),cc2h(:,2));
    rho_m = rho_m*mass/sum(rho_m.*area2h);
    rhoa = RHt*It*rho_all;
    [rhos]=compute_rhosigma(ind,edges,cc,mid,N,rho_f,rho_in,gradt,Mst,RHt,It,Rst,rec,uk,'rhos');
    
    [~,indc]=min((cc(:,1)-0.5).^2+(cc(:,2)-0.5).^2);
    phi = phi+(potential(0.5,0.5,0)-phi(indc));
    
    % Compute the errors
    err_cost(h_i) = abs(W2-W2th); % cost error
    err_mid(h_i) = sum(Mx2h*abs(rhomid-rho_m)); % density error in the midpoint
    for k=1:Nt
        t = (k-0.5)/Nt;
        phi_real = potential(cc(:,1),cc(:,2),t);
        gphi_real = grad*phi_real;
        rho_real = geodesic(cc2h(:,1),cc2h(:,2),t);
        rho_real = rho_real*mass/sum(rho_real.*area2h);
        rho_t = 0.5*(rho_all((k-1)*ncell2h+1:k*ncell2h)+rho_all(k*ncell2h+1:(k+1)*ncell2h));
        % error on the potential
        err_p(h_i) = err_p(h_i) + (1/Nt)*sum( rhoa((k-1)*ncell+1:k*ncell).*Mx*(phi_real-phi((k-1)*ncell+1:k*ncell)).^2 );
        % error on the gradient of the potential
        err_gp(h_i) = err_gp(h_i) + (1/Nt)*sum( rhos((k-1)*nei+1:k*nei).*Ms*(gphi_real-gradphi((k-1)*nei+1:k*nei)).^2);
        % error on the geodesic
        err_rhoth(h_i) = err_rhoth(h_i) + (1/Nt)*sum(Mx2h*abs(rho_t-rho_real));
    end
    err_p(h_i) = sqrt(err_p(h_i));
    err_gp(h_i) = sqrt(err_gp(h_i));
    
    
    % print mid point
    fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rhomid,'FaceColor','flat');
    colorbar
    axis equal
    axis off
    str =num2str(h_i);
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


