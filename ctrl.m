clear
close all


% Computation of Optimal Transport with finite volumes
% (https://hal.archives-ouvertes.fr/hal-03032446)

% Time discretization. Uniform grid with Nt=N+1 intervals. N is the number
% of intermediate densities. Choose N odd, N>=1, in order to have the
% approximate geodesic midpoint.
N = 100;
Nt = N+1;

% Space discretization. Three types of mesh families available:
% 1 -> regular triangulation of the domain, with only acute angles
%      (https://www.i2m.univ-amu.fr/fvca5/benchmark/Meshes/index.html)
% 2 -> nested-constructed meshes, based on the previous family
%      (find details at https://hal.archives-ouvertes.fr/hal-03032446)
% 3 -> cartesian grids
% For each mesh, five levels of refinement h_i, 1->5, are available.
mesh_type = 2;
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


grounded_node=1
solver_approach=8
plot=0

ctrl_inner=ctrl_solver;
ctrl_outer=ctrl_solver;
ctrl_inner.init('\',1e-5,3000,0.1);
ctrl_outer.init('jacobi',1e-5,3000);
compute_eigen=0
