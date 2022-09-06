# Computation of Optimal Transport with finite volumes. 

Approximate solution of the dynamical form of the L2 optimal trasport problem with finite volumes discretization. This code implements the scheme presented in:

[1] A. Natale, G. Todeschi, "Computation of Optimal Transport with Finite Volumes", ESAIM: Mathematical Modelling and Numerical Analysis, 55(5):1847-1871, 2021, https://www.esaim-m2an.org/articles/m2an/abs/2021/06/m2an210008/m2an210008.html.

and the work on the linear algrebra parts described in:

[2] E. Facca, G. Todeschi, A. Natale, M. Benzi "Efficient preconditioners for solving dynamical optimal transport via interior point methods". Preprint available arXiv https://arxiv.org/abs/2209.00315.

The file "l2otp_solve.m" contains the main solver taking the meshes, the initial and the final densities ot be transported, and the solver controls and computing the optimal interpolationg density and the optimal potential.
Input and output data are enclosed in adhoc structures.
Examples of usage is contained in the "run_NataleTodeschi.m"  reproducing the experiments in [1] and 
"run_FaccaTodeschiNataleBenzi.m" file, which reproduces the linear algebra experiments described in [2].
Note that, in order to fully exploit the linear algebra capabilities of the solver, you need to install the multigrid solver [AGMG](http://agmg.eu/) solver (see agmg/README.md first).

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

