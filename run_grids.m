twogrids = 1;
fn1 = 'meshes/mesh_gen/tri2/mesh1/coord.txt'; %coordinates
fn2 = 'meshes/mesh_gen/tri2/mesh1/topol.txt'; %topology
[nodes,cells] = get_mesh(fn1,fn2,0,0,0,0);
grid_rho = TPFA_grid;
grid_rho.init(cells,nodes);
if twogrids
	%grid_phi = TPFA_grid;
	[grid_phi, I] = grid_rho.refine();
else
	grid_phi = grid_rho;
	I = speye(grid_rho.ncell);
end

[grid_rho_right, grid_phi_right, I_right] = init_grids(2, 1); % 2=not reordered, 1= coarse mesh


fprintf('difference edges grid_rho %f\n', norm(grid_rho.edges(:,1:7)-grid_rho_right.edges(:,1:7)))
%[grid_rho.cells,grid_rho_right.cells]

fprintf('difference ind grid_rho %f\n', norm(grid_rho.ind.internal-grid_rho_right.ind.internal))

disp('grid_phi')
fprintf(' ncell : new %d - old %d \n',grid_phi.ncell,grid_phi_right.ncell)
fprintf(' max index cells : new %d - old %d | right %d \n',max(max(grid_phi.cells(:,2:5))),max(max(grid_phi_right.cells(:,2:5))),grid_rho.ncell+grid_rho.nedges+grid_rho.nnodes)
fprintf('difference cells grid_rho %f\n', norm(grid_phi.cells(:,1:5)-grid_phi_right.cells(:,1:5)))
%[grid_phi.cells,grid_phi_right.cells]
fprintf(' nnode : new %d - old %d \n',grid_phi.nnodes,grid_phi_right.nnodes)

fprintf(' nedge : new %d - old %d \n',grid_phi.nedges,size(grid_phi_right.edges,1))
fprintf(' nedge internal : new %d old %d \n',grid_phi.nei,grid_phi_right.nei)
