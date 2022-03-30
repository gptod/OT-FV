% class spatial information of grid
classdef TPFA_grid <handle
  properties
		ncell;
    nodes;
    cells;
    edges;
    ind;
    area;
    cc;
    mid;
    h;
  end
  methods
		% constructor
    function obj = init(cells,nodes)
		  obj.cells=cells;
			obj.node=nodes;
    end
   	% destructor
		function obj = kill(obj)
			if (obj.initalized == 1)
				clear obj.cells;
				clear obj.edges
				clear obj.ind;
				clear obj.area;
				clear obj.cc;
				clear obj.mid;
			end
		end

		function [mass,edge_mass,div,grad] = build_matrices(obj);
			mass = spdiags(obj.area,0,obj.ncell,obj.ncell);
			nei=size(obj.ind.internal,1);
			ds = obj.edges(obj.ind.internal,5).*obj.edges(obj.ind.internal,6);
			edge_mass = spdiags(ds,0,nei,nei);
			div = Div2D(obj.ncell,nei,obj.ind,obj.edges); % divergence matrix
			grad = -edge_mass\div'; % gradient matrix
		end
		
		% info
		function obj = info(obj,fid)
      if (~exist('fid','var') )
				fid=1;
			end
		end
	end
end
