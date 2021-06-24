clear
close all


% Convergence test
% test 1: sin -> pure translation
% test 2: compression -> contraction
test_case='cross';

% Type of mesh
mesh_type = 2;

% Type of reconstruction
% 1 -> weighted arithmetic mean
% 2 -> weighted harmonic mean
rec = 2;


eps_0 = 1e-6; % tolerance 

verb = 0;

plot = 0;
plot_figures=0;
restart=0;
save_data=0;
compute_eigen=0;

verbose=0;



h_i = 1;
N = 2*2^(h_i-1)-1;
%N = 4*2^(h_i-1)-1;
Nt = N+1;

str_test=sprintf('%s_h%d_rec%d_N%0.5d_',test_case,h_i,rec,N)


solver_approach=11;

ctrl_inner11=ctrl_solver;
ctrl_inner22=ctrl_solver;
ctrl_outer=ctrl_solver;

ctrl_outer.init('fgmres',1e-5,3000,0.0,0); % verbose=[1,2] works only for fgmres
%ctrl_outer.init('fgmres',1e-5,3000,0.0,0);
left_right='right';

outer_prec='full';

% set here other approximate inverse of block11
ctrl_inner11.init('diag',... %approach
        1e-12,... %tolerance
        10,...% itermax
        0.0,... %omega
        0,... %verbose
        'diag'); %label
extra_info='full';
%extra_info='block';
relax4_inv11=1e-12;

% set grounded_node>0 to gorund the potential in grounded node
grounded_node=0;

relax4_inv22=0;

ctrl_inner22.init('agmg',1e-13,10,1.0,0,'agmg10');
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


save('save/rho_cross_subtri2_harm_h5N32','rho')


if plot==1
figure
fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rho_in,'FaceColor','flat');
%colorbar
axis square
axis off
caxis([0 max(rho_all)])
str =num2str(0);
outname=strcat('figures/cross_subtri2_harm_h5N32_rhok',str,'.jpg');
saveas(fig,outname)
figure
fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rho_f,'FaceColor','flat');
%colorbar
axis square
axis off
caxis([0 max(rho_all)])
str =num2str(N+1);
outname=strcat('figures/cross_subtri2_harm_h5N32_rhok',str,'.jpg');
saveas(fig,outname)
dsp = 1;
%dsp = ceil(N/2);
for k=dsp:dsp:N-dsp+1
    rhok = rho((k-1)*ncell2h+1:k*ncell2h);
    figure
    fig=patch('Faces',cells2h(:,2:end),'Vertices',nodes2h,'edgecolor','none','FaceVertexCData',rhok,'FaceColor','flat');
    %colorbar
    axis square
    axis off
    caxis([0 max(rho_all)])
    str =num2str(k);
    outname=strcat('figures/cross_subtri2_harm_h5N32_rhok',str,'.jpg');
    saveas(fig,outname)
end
end