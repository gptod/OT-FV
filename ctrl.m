%clear
%close all


% Computation of Optimal Transport with finite volumes
% (https://hal.archives-ouvertes.fr/hal-03032446)

% Time discretization. Uniform grid with Nt=N+1 intervals. N is the number
% of intermediate densities. Choose N odd, N>=1, in order to have the
% approximate geodesic midpoint.
%N = 16;
Nt = N+1;

%test_case='gauss'
test_case='sin'
%test_case='sin_old'
%test_case='cross'


% Space discretization. Three types of mesh families available:
% 1 -> regular triangulation of the domain, with only acute angles
%      (https://www.i2m.univ-amu.fr/fvca5/benchmark/Meshes/index.html)
% 2 -> nested-constructed meshes, based on the previous family
%      (find details at https://hal.archives-ouvertes.fr/hal-03032446)
% 3 -> cartesian grids
% For each mesh, five levels of refinement h_i, 1->5, are available.
mesh_type = 2;
h_i = 2;
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


grounded_node=1
% 9 : P=(A B1T)
%       (0 -C )

% 10 : P=(~S B1T)
%        (0  -C )
% ~S=A+B1T C^{-1} B2 
% ~S=upper_triang( A+B1T C^{-1} B2 )
% ~S=block_diag(   A+B1T C^{-1} B2 ) 

solver_approach=11
			
plot=0

ctrl_inner11=ctrl_solver;
ctrl_inner22=ctrl_solver;
ctrl_outer=ctrl_solver;


%ctrl_inner11.init('krylov',1e-07,1000,1.0,1);
%ctrl_inner11.init('agmg',1e-12,1000,1.0,0);
%ctrl_outer.init('1',1e-14,3000);
ctrl_inner11.init('diag',1e-12,1,1.0,0);
extra_info='block_triang'

ctrl_inner22.init('incomplete',1e-12,10,1.0,0);


ctrl_outer.init('fgmres',1e-6,3000);


compute_eigen=0;
verbose=1;


str_test=sprintf('%s_h%d_rec%d_N%0.5d_',test_case,h_i,rec,N)
disp(str_test)

controls = struct('indc',grounded_node,...
		  'sol',solver_approach,...
		  'ctrl_inner11',ctrl_inner11,...
		  'ctrl_inner22',ctrl_inner22,...
		  'ctrl_outer',ctrl_outer,...
		  'compute_eigen',compute_eigen,...
		  'verbose',verbose,...
		  'extra_info',extra_info)

approach_string=def_string_approach(controls)
