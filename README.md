# Computation of Optimal Transport with finite volumes. 
prova
Approximate solution of the dynamical form of the L2 optimal trasport problem with finite volumes discretization.
This code implements the scheme presented in: A. Natale, and G. Todeschi, "Computation of Optimal 
Transport with Finite Volumes", arxiv preprint, arXiv:2012.00349 (https://arxiv.org/abs/2012.00349).

The files "geod.m" and "convergence.m" are the main sources. "geod.m" is the code 
which can be used to perform general tests (as the ones contained in the paper),
whereas "convergence.m" is the specific code which performs the convergence tests
of the paper. In the former case, the function "bc_density" can be used to change
the boundary conditions on the density, i.e. the initial and final density. 
In the latter case, only the two tests presented in the paper are available,
but other can be equivalently added (pay attention at the global constant for
the potential though!).

The folder "meshes" contains the meshes data structures used by the two aforementioned codes.
Three types of mesh families are available: \
1 -> regular triangulation of the domain, with only acute angles
     (https://www.i2m.univ-amu.fr/fvca5/benchmark/Meshes/index.html) \
2 -> nested-constructed meshes, based on the previous family
     (find details in the paper) \
3 -> cartesian grids \
For each mesh, five levels of refinement h_i, 1->5, are available. \
Mesh structure: \
nodes -> array of nodes coordinates [x y] \
cells -> array of cells nodes [#nodes node1 node2 node3 ...] \
edges -> array of edges [node1 node2 K L d_sigma m_sigma m_sigma/d_sigma] \
ind -> indices struct: ind.internal indices of internal edges, ind.bound indices of boundary edges \
cc -> array of cell centers coordinates [x y] \
area -> array of cell measures \
mid -> array of midpoints coordinates [x y] \
h -> meshsize, max(diam(K)) or max(sqrt(area)) \
The outer mesh structure, if any, is labelled with 2h to distinguish it
from the inner nested one.

The mesh type, the level of space and time refinement, the type of reconstruction
used for the finite volume scheme, can be selected via the specific parameters 
available at the beginning of the source codes.

