function [grid_rho,grid_phi,I]=init_grids_3d(mesh_type,h_i)

	%% read mesh
	if mesh_type == 1
        mesh_name = strcat('meshes/cube',num2str(h_i));
	else
		disp('Mesh not supprted')
		return
  end


	% load mesh
	clear 'I'
	mesh=load(mesh_name);
	

	
	% init objects and assign "manually" the class member
	grid_phi = TPFA_grid;
	grid_phi.ncells = mesh.ncells;
	grid_phi.nodes = mesh.nodes;
	grid_phi.nnodes = size(mesh.nodes,1);
	grid_phi.cells = mesh.cells;
  grid_phi.ncells = mesh.ncells;
	grid_phi.edges = mesh.edges;
  grid_phi.nedges = mesh.nedges;
  grid_phi.faces = mesh.faces;
  grid_phi.nfaces = mesh.nfaces;
  grid_phi.sigma = mesh.sigma;
	grid_phi.ind   = mesh.ind;
	grid_phi.area  = mesh.area;
	grid_phi.cc    = mesh.cc;
  grid_phi.mid_faces   = mesh.mid_faces;
  grid_phi.nsig = mesh.nsig;
	grid_phi.nsig_in = size(mesh.ind.internal,1);

	ncells=grid_phi.ncells;


	% if the injection operator I does not exist in the mesh structure,
	% the nested sub-mesh coincides with the outer one (labeled with 2h)
	
	if isfield(mesh,'I')
		grid_rho = TPFA_grid;

		% assign "manually" the class member
		grid_rho.ncells = mesh.ncells2h;
		grid_rho.nodes = mesh.nodes2h;
    grid_rho.nnodes = mesh.nnodes2h;
		grid_rho.cells = mesh.cells2h;
		grid_rho.edges = mesh.edges2h;
    grid_rho.nedges = mesh.nedges2h;
    grid_rho.faces = mesh.faces2h;
    grid_rho.nfaces = mesh.nfaces2h;
    grid_rho.sigma = mesh.sigma2h;
		grid_rho.ind   = mesh.ind2h;
		grid_rho.area  = mesh.area2h;
		grid_rho.cc    = mesh.cc2h;
		grid_rho.mid_faces   = mesh.mid_faces2h;
    grid_rho.nsig = mesh.nsig2h;
		grid_rho.nsig_in = mesh.nsig_in2h;
		
		I=mesh.I;
	else
		I = speye(ncells,ncells);
		
		% init object just by copy
		grid_rho = grid_phi;

	end

end
