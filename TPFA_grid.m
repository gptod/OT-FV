% class spatial information of grid
classdef TPFA_grid <handle
  properties
		ncell;
		cells;
		nodes;
    nnodes;
    edges;
    ind;
    area;
    cc;
    mid;
    h;
  end
  methods
		% constructor from topology and coordinates
		% INPUTS
		% topol : interger (nnodeincell, ncell) array with mesh topology
		%         nnodeincell = 3, 4
		% coord : real(2,nnode) array with nodes coordinates
		% OUTPUTS
		% obj: mesh initialized
    function obj = init(topol,coord)
			% add the number of node in cells
			% TODO: use csr matrix
		  obj.cells = [sum(topol~=0,2) topol];
			obj.ncell = size(obj.cells,2);
			obj.node = nodes;
			obj.nnodes = size(nodes,2);

			% compute cells related properties
			[obj.cc,obj.area,h] = mesh_properties(obj.nodes,obj.cells);

			% compute edges structure
			[obj.edges,obj.mid] = str_edge(obj.nodes, obj.cells, obj.cc);
			nedge = size(edges,1);

			% compute indices 
			obj.ind = indices(obj.ncell,obj.edges);
    end

		% Generate a finer mesh subdividiving each cell
		% in 3 (triangular mesh) or 4 subcells (quadrilater cells).
		% 
		% INPUTS:
		% obj : mesh (initialized)
		% OUPUTS:
		% refined : refined mesh
		% coarse2fine : sparse matrix mapping functions/vectors on the coarse
		%               mesh into functions/vectors in the fine one.
		%               It also describes the relation between cooarse and fine cells
		function refined, coarse2fine = sub_mesh(obj)
			% shorthands
			ncell = obj.ncell;
			nedge = obj.nedge;
			nnode = obj.nnode;
			

			
			% compute cells-edges structure
			maxn= max(obj.cells(:,1));
			[cell_e,~,~] = str_cell(maxn,obj.ind,obj.edges,obj.mid,obj.cc);
		
			cells_f = zeros(3*ncell,5);
			nodes_f = zeros(nnode+nedge+ncell,2);
			% TODO: remove this, derived variables are computed in the constructor
			nodes = [obj.nodes zeros(nnode,1)];
			mid = [obj.mid zeros(nedge,1)];
			cc = [obj.cc zeros(ncell,1)];
			int = 0;
			coarse2fine = sparse(3*ncell,ncell);

			for k=1:ncell
				
				K = obj.cells(k,2:4);
				E = cell_e(k,2:4);
				
				for j=1:3
					
					
					ej = E(obj.edges(E',1)==K(j) | obj.edges(E',2)==K(j));
					
					if obj.nodes(K(j),3)==0
            int = int+1;
            A = int;
            obj.nodes(K(j),3) = int;
            nodes_f(int,:) = obj.nodes(K(j),1:2);
					else
            A = obj.nodes(K(j),3);
					end
					
					if mid(ej(1),3)==0
            int = int+1;
            B = int;
            mid(ej(1),3) = int;
            nodes_f(int,:) = mid(ej(1),1:2);
					else
            B = mid(ej(1),3);
					end
					
					if cc(k,3)==0
            int = int+1;
            C = int;
            cc(k,3) = int;
            nodes_f(int,:) = cc(k,1:2);
					else
            C = cc(k,3);
					end
					
					if mid(ej(2),3)==0
            int = int+1;
            D = int;
            mid(ej(2),3) = int;
            nodes_f(int,:) = mid(ej(2),1:2);
					else
            D = mid(ej(2),3);
					end
					
					cells_f(3*(k-1)+j,:) = [4 A B C D];
					
				end
				
				coarse2fine(3*(k-1)+[1;2;3],k) = ones(3,1);
				
			end

			refined = TPFA_grid;
			refined.init(cells_f,coord_f);
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


		% create the matrices used in the finite volume discretization
		% INPUTS
		% obj : mesh
		% OUTPUTS:
		% mass : diagonal matrix with cell volume
		% edge_mass : diagonal matrix with edge length
		% div : sparse matrix discretizing the divergence operator
		% grad : sparse matrix discretizing the grad operator
		%        (grad = - edge_mass^{-1} div )
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

		function [cc,area,h] = mesh_properties(nodes,cells)

			% INPUT:
			% node: nodes' coordinate
			% cell: topology
			% OUTPUT:
			% cc: circumcenters of the polygons
			% area: area of the polygons
			% h: meshsize of the mesh


			% compute the circumcenters of the cells, namely the intersection of the
			% segments orthogonal to the edges and passing through their midpoints
			% take for this purpose two random edges...
			vert = cells(:,2:end); % disregard the column specifying the number of vertices per cell
      % attention: the vert mat contains possibly something
      % different from vertices at the right end!!!
			nvert = cells(:,1);

			% compute the orthogonal vector to the first edge vertex_1-vertex_2
			v = [nodes(vert(:,1),2)-nodes(vert(:,2),2) nodes(vert(:,2),1)-nodes(vert(:,1),1)];
			% use the second edge vertex_2-vertex_3 if a node's x-ccordinate is equal to zero
			corrv = find(v(:,1)==0);
			v(corrv,:) = [nodes(vert(corrv,2),2)-nodes(vert(corrv,3),2) nodes(vert(corrv,3),1)-nodes(vert(corrv,2),1)];
			if any(v(:,1)==0)
				error('there is a problem with the mesh')
			end

			% compute the orthogonal vector to the last edge vertex_1-vertex_nvert
			w = zeros(size(vert,1),2);
			for i=1:size(vert,1)
				w(i,:) = [nodes(vert(i,1),2)-nodes(vert(i,nvert(i)),2) nodes(vert(i,nvert(i)),1)-nodes(vert(i,1),1)];
			end
			% detect parallel edges (v//w) and change the edge w: use instead the edge
			% vertex_(nvert-1)-vertex_nvert 
			corrw = [];
			for i=1:size(v)
				if w(i,:)-v(i,:)==[0 0]
					corrw = [corrw; i];
				end
			end
			for i=1:length(corrw)
				w(corrw(i),:) = [nodes(vert(corrw(i),nvert(corrw(i))-1),2)-nodes(vert(corrw(i),nvert(corrw(i))),2) nodes(vert(corrw(i),nvert(corrw(i))),1)-nodes(vert(corrw(i),nvert(corrw(i))-1),1)];
			end

			% compute the midpoints of the two edges considered
			midv = [0.5*(nodes(vert(:,1),1)+nodes(vert(:,2),1)) 0.5*(nodes(vert(:,2),2)+nodes(vert(:,1),2))];
			midv(corrv,:) = [0.5*(nodes(vert(corrv,2),1)+nodes(vert(corrv,3),1)) 0.5*(nodes(vert(corrv,3),2)+nodes(vert(corrv,2),2))];
			midw = zeros(size(vert,1),2);
			for i=1:size(vert,1)
				midw(i,:) = [0.5*(nodes(vert(i,1),1)+nodes(vert(i,nvert(i)),1)) 0.5*(nodes(vert(i,nvert(i)),2)+nodes(vert(i,1),2))];
			end
			for i=1:length(corrw)
				midw(corrw(i),:) = [0.5*(nodes(vert(corrw(i),nvert(corrw(i))-1),1)+nodes(vert(corrw(i),nvert(corrw(i))),1)) 0.5*(nodes(vert(corrw(i),nvert(corrw(i))),2)+nodes(vert(corrw(i),nvert(corrw(i))-1),2))];
			end

			% compute the point of intersection between the two vectors
			beta = (midw(:,2)-midv(:,2)-v(:,2).*(midw(:,1)-midv(:,1))./v(:,1))./(v(:,2).*w(:,1)./v(:,1)-w(:,2));
			cc = midw + beta.*w;


			% a = [nodes(cells(:,1),1) nodes(cells(:,1),2) nodes(cells(:,1),3)];
			% b = [nodes(cells(:,2),1) nodes(cells(:,2),2) nodes(cells(:,2),3)];
			% c = [nodes(cells(:,3),1) nodes(cells(:,3),2) nodes(cells(:,3),3)];
			% ac = c-a;
			% ab = b-a;
			% abxac = cross(ab,ac,2);
			% abxacxab = cross(abxac,ab,2);
			% acxabxac = cross(ac,abxac,2);
			% cc = zeros(size(tri,1),3);
			% for k=1:size(tri,1)
			%     cc(k,:) = a(k,:) + (norm(ac(k,:),2)^2*abxacxab(k,:) + norm(ab(k,:),2)^2*acxabxac(k,:))/(2*norm(abxac(k,:),2)^2);
			% end


			% compute the area of the polygons
			area = zeros(size(vert,1),1);
			for i=1:size(vert,1)
				vert_i = vert(i,1:nvert(i));
				x_i = nodes(vert_i,:);
				area(i) = area_pol(nvert(i),x_i);
			end

			% for a manifold:
			% for k=1:size(tri,1)
			%     x1 = node(tri(k,1),1);
			%     x2 = node(tri(k,2),1);
			%     x3 = node(tri(k,3),1);
			%     y1 = node(tri(k,1),2);
			%     y2 = node(tri(k,2),2);
			%     y3 = node(tri(k,3),2);
			%     z1 = node(tri(k,1),3);
			%     z2 = node(tri(k,2),3);
			%     z3 = node(tri(k,3),3);
			%     areak = 0.5*norm(cross([x2-x1; y2-y1; z2-z1],[x3-x1; y3-y1; z3-z1]),2);
			%     area = [area; areak];
			% end

			% compute the meshsize: max(diam(polygon))
			if any(nvert~=3)
				% to simplify the problem use the square root of the area instead of
				% the diameter 
				%h = max(sqrt(area));
				
				% for a convex polygon the diameter is equal to the maximum distance
				% between vertices
				h = 0;
				for i=1:size(vert,1)
					for j=1:nvert(i)
            for k=j+1:nvert(i)
              h = max(h, sqrt((nodes(vert(i,j),1)-nodes(vert(i,k),1)).^2+(nodes(vert(i,j),2)-nodes(vert(i,k),2)).^2) );
            end
					end
				end
			else
				% if the polygons are triangle
				h1 = sqrt((nodes(vert(:,1),1)-nodes(vert(:,2),1)).^2+(nodes(vert(:,1),2)-nodes(vert(:,2),2)).^2);
				h2 = sqrt((nodes(vert(:,1),1)-nodes(vert(:,3),1)).^2+(nodes(vert(:,1),2)-nodes(vert(:,3),2)).^2);
				h3 = sqrt((nodes(vert(:,2),1)-nodes(vert(:,3),1)).^2+(nodes(vert(:,2),2)-nodes(vert(:,3),2)).^2);
				h = max(max([h1 h2 h3]));
			end

		end

		function [edges,mid] = str_edge(nodes,cells,cc)

			% rough computation of the edges structure
			% INPUT:
			% node: nodes' coordinates
			% cell: topology
			% cc: cell orthocenters
			% OUTPUT:
			% edges: edges structure
			% mid: midpoints vector

			vert = cells(:,2:end); % disregard the column specifying the number of vertices per cell

			nvert = cells(:,1);
			nnode = size(nodes,1);
			
			% construct all the edges with repetitions
			% the edge is construct as: [vertexA vertexB], vertexA and vertexB are vertices numbers,
			% with vertexA < vertexB
			edges_rep = zeros(sum(nvert),3);
			k = 0;
			for i=1:size(vert,1)
				for j=1:nvert(i)
					k = k+1;
					if j==nvert(i)
            edges_rep(k,:) = [vert(i,j) vert(i,1) i];
					else
            edges_rep(k,:) = [vert(i,j) vert(i,j+1) i];
					end
				end
			end
			
			% order each edge from the smaller vertex to the bigger one
			[edges_rep(:,1:2),~] = sort(edges_rep(:,1:2),2);
			% order the edges in ascending order with respect to their first vertex
			[~,ord2] = sort(edges_rep(:,1),1);
			edges_rep = edges_rep(ord2,:);

			% re-order the edges in ascending order also with respect to their second vertex
			row = 1;
			for i=1:nnode
				k=1;
				if row+k+1>size(edges_rep)
					break
				end
				while edges_rep(row+k,1)==edges_rep(row,1)
					k = k+1;
					if row+k>size(edges_rep)
            break
					end
				end
				[~,ord2]=sort(edges_rep(row:row+k-1,2),1);
				A = edges_rep(row:row+k-1,:);
				edges_rep(row:row+k-1,:) = A(ord2,:);
				row = row+k;
			end

			% construct the edges structure without repetitions
			edges = [edges_rep(1,:) 0];
			row = 1;
			for i=2:size(edges_rep,1)
				if edges_rep(i,1:2)==edges(row,1:2)
					edges(row,4) = edges_rep(i,3);
					%edges(row,5) = edges(row,5)+edges_rep(i+1,5);
				else
					edges = [edges; [edges_rep(i,:) 0]];
					row = row+1;
				end
			end

			% compute edges midpoints
			mid = [0.5*(nodes(edges(:,1),1)+nodes(edges(:,2),1)) 0.5*(nodes(edges(:,1),2)+nodes(edges(:,2),2))];


			% add the edges' distances d_sig, length m_sig and the transmissivity d_sig/m_sig
			% add also the distances d_{K,sig}, d_{L,sig}
			edges = [edges zeros(size(edges,1),5)];
			for e=1:size(edges,1)
				if edges(e,4)~=0
					d_sig = sqrt( (cc(edges(e,3),1)-cc(edges(e,4),1)).^2 + (cc(edges(e,3),2)-cc(edges(e,4),2)).^2 );
					dKsig = sqrt((mid(e,1)-cc(edges(e,3),1)).^2+...
											 (mid(e,2)-cc(edges(e,3),2)).^2);
					dLsig = sqrt((mid(e,1)-cc(edges(e,4),1)).^2+...
											 (mid(e,2)-cc(edges(e,4),2)).^2);
				else
					d_sig = sqrt( (cc(edges(e,3),1)-mid(e,1)).^2 + (cc(edges(e,3),2)-mid(e,2)).^2 );
					dKsig = sqrt((mid(e,1)-cc(edges(e,3),1)).^2+...
											 (mid(e,2)-cc(edges(e,3),2)).^2);
					dLsig = 0;
				end
				m_sig = sqrt((nodes(edges(e,1),1)-nodes(edges(e,2),1)).^2+(nodes(edges(e,1),2)-nodes(edges(e,2),2)).^2);
				edges(e,5:9) = [d_sig m_sig m_sig/d_sig dKsig dLsig];
			end



		end

		function ind = indices(ncell,edges)

			% INPUT:
			% ncell: number of cells
			% edges: edges structure
			% OUPUT:
			% ind: indices structure

			% this function creates indices for a simplified access to the matrices
			% these indices enable to select the entries of submatrices of the FV
			% matrices

			% indices for internal edges:

			internal = find(edges(:,4)~=0);
			% ind.i_KK = sub2ind([ncell ncell],edges(internal,3),edges(internal,3));
			% ind.i_LL = sub2ind([ncell ncell],edges(internal,4),edges(internal,4));
			% ind.i_KL = sub2ind([ncell ncell],edges(internal,3),edges(internal,4));
			% ind.i_LK = sub2ind([ncell ncell],edges(internal,4),edges(internal,3));

			% [ind.u_i,~,ind.o_i] = unique(edges(internal,3));
			% [ind.u_e,~,ind.o_e] = unique(edges(internal,4));

			% [ind.u_KK,~,ind.o_KK] = unique(ind.i_KK);
			% [ind.u_LL,~,ind.o_LL] = unique(ind.i_LL);
			% [ind.u_KL,~,ind.o_KL] = unique(ind.i_KL);
			% [ind.u_LK,~,ind.o_LK] = unique(ind.i_LK);

			ind.internal = internal;
			ind.bound = find(edges(:,4)==0);

			ind.nei = length(internal);
			% ind.et_K = sub2ind([ind.nei ncell],[1:ind.nei]',edges(internal,3));
			% ind.et_L = sub2ind([ind.nei ncell],[1:ind.nei]',edges(internal,4));

			%[ind.etu_K,~,ind.eto_K] = unique(ind.et_K);
			%[ind.etu_L,~,ind.eto_L] = unique(ind.et_L);

			% indices for all edges:

			% [ind.allu_i,~,ind.allo_i] = unique(edges(:,3));
			% [ind.allu_e,~,ind.allo_e] = unique(edges(internal,4));
			% 
			% ind.alli_ii = sub2ind([ncell ncell],edges(:,3),edges(:,3));
			% ind.alli_ee = sub2ind([ncell ncell],edges(internal,4),edges(internal,4));
			% ind.alli_ie = sub2ind([ncell ncell],edges(internal,3),edges(internal,4));
			% ind.alli_ei = sub2ind([ncell ncell],edges(internal,4),edges(internal,3));
			% 
			% [ind.allu_ii,~,ind.allo_ii] = unique(ind.alli_ii);
			% [ind.allu_ee,~,ind.allo_ee] = unique(ind.alli_ee);
			% [ind.allu_ie,~,ind.allo_ie] = unique(ind.alli_ie);
			% [ind.allu_ei,~,ind.allo_ei] = unique(ind.alli_ei);
			% 
			% nedge = size(edges,1);
			% ind.allet_i = sub2ind([nedge ncell],[1:nedge]',edges(:,3));
			% ind.allet_e = sub2ind([nedge ncell],internal,edges(internal,4));
			% [ind.alletu_i,~,ind.alleto_i] = unique(ind.allet_i);
			% [ind.alletu_e,~,ind.alleto_e] = unique(ind.allet_e);

		end

		function [cell_e,cell_eint,cell_dist] = str_cell(maxn,ind,edges,mid,cc)

			% INPUT:
			% maxn: maximum number of vertices/edges per polygonal cell
			% ind: indices structure
			% edges: edges structure
			% mid: edges midpoints
			% cc: othocenters of the cells
			% OUTPUT:
			% cell_e: edges per cell structure
			% cell_eint: internal edges per cell structure
			% cell_dist: distances from cc to the edges per cell structure

			ncell = size(cc,1);
			nedge = size(edges,1);
			nei = length(ind.internal);

			cell_eint = zeros(ncell,maxn+1);
			% for each internal edge e, take the two cells K and L defining it and assign to both the
			% edge e
			for e=1:nei
				K = edges(ind.internal(e),3);
				L = edges(ind.internal(e),4);
				cell_eint(K,1) = cell_eint(K,1)+1;
				cell_eint(K,cell_eint(K,1)+1) = e;
				if L~=0
					cell_eint(L,1) = cell_eint(L,1)+1;
					cell_eint(L,cell_eint(L,1)+1) = e;
				end
			end

			cell_e = zeros(ncell,maxn+1);
			% for each edge e, take the two cells K and L defining it and assign to both the
			% edge e
			for e=1:nedge
				K = edges(e,3);
				L = edges(e,4);
				cell_e(K,1) = cell_e(K,1)+1;
				cell_e(K,cell_e(K,1)+1) = e;
				if L~=0
					cell_e(L,1) = cell_e(L,1)+1;
					cell_e(L,cell_e(L,1)+1) = e;
				end
			end


			cell_dist = zeros(ncell,maxn+1);
			cell_dist(:,1) = cell_e(:,1);
			% for each cell i, for each edges j of the cell i, compute the distance |cc(i)-mid(j)|
			for i=1:ncell
				for j=1:cell_dist(i,1)
					cell_dist(i,1+j) = sqrt((cc(i,1)-mid(cell_e(i,1+j),1)).^2+(cc(i,2)-mid(cell_e(i,1+j),2)).^2);
				end
			end

		end
		
	end
end
