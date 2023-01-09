function [grid_rho,grid_phi,I]=init_grids(mesh_type,h_i)

	%% read mesh
	if mesh_type == 1
    mesh_name = strcat('meshes/tri2_mesh',num2str(h_i));
	elseif mesh_type == 2
    mesh_name = strcat('meshes/subtri2_mesh',num2str(h_i));
	elseif mesh_type == 3
		mesh_name = strcat('meshes/sq_mesh',num2str(h_i));
	elseif mesh_type == 4
    mesh_name = strcat('meshes/tri2ord_mesh',num2str(h_i));
	elseif mesh_type == 5
		mesh_name = strcat('meshes/subtri2ord_mesh',num2str(h_i));
	elseif mesh_type == 6
		mesh_name = strcat('meshes/sqord_mesh',num2str(h_i));
	else
		disp('Mesh not supprted')
		return
	end


	% load mesh
	clear 'I'
	mesh=load(mesh_name);
	

	
	% init objects and  assign "manually" the class member
	grid_phi = TPFA_grid;
	grid_phi.ncell = mesh.ncell;
	grid_phi.nodes = mesh.nodes;
	grid_phi.cells = mesh.cells;
	grid_phi.edges = mesh.edges;
	grid_phi.ind   = mesh.ind;
	grid_phi.area  = mesh.area;
	grid_phi.cc    = mesh.cc;
	grid_phi.mid   = mesh.mid;
	grid_phi.nei = size(mesh.ind.internal,1);

	ncell=grid_phi.ncell;


	% if the injection operator I does not exist in the mesh structure,
	% the nested sub-mesh coincides with the outer one (labeled with 2h)
	
	if isfield(mesh,'I')
		grid_rho = TPFA_grid;

		% assign "manually" the class member
		grid_rho.ncell = mesh.ncell2h;
		grid_rho.nodes = mesh.nodes2h;
		grid_rho.cells = mesh.cells2h;
		grid_rho.edges = mesh.edges2h;
		grid_rho.ind   = mesh.ind2h;
		grid_rho.area  = mesh.area2h;
		grid_rho.cc    = mesh.cc2h;
		grid_rho.mid   = mesh.mid2h;

		grid_rho.nei = size(mesh.ind2h.internal,1);
		
		I=mesh.I;
	else
		I = speye(ncell,ncell);
				
		% init object just by copy
		grid_rho = grid_phi;

	end

end
