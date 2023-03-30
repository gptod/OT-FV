% class spatial information of grid
classdef TPFA_grid <handle
  properties
	% cells, cells number
	ncells;
	cells;
	% nodes coordinate, nodes number
	nodes;
    nnodes;
	% edges, edge, edge number, internal edge number
    edges;
	nedges;
    faces;
    nfaces;
    sigma;
    nsig;
	nsig_in;
    nsig_b;
    cell_sig;
    ind;
    area;
    cc;
    mid_edges;
    mid_faces;
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
    function obj = init(obj, topol,nodes)
		% add the number of node in cells
		% TODO: use csr matrix

	    obj.cells = [sum(topol~=0,2) topol];
		obj.ncells = size(obj.cells,1);
		obj.nodes = nodes;
		obj.nnodes = size(nodes,1);

		% compute cells related properties
		[obj.cc,obj.area,obj.h] = obj.mesh_2d(obj.nodes,obj.cells);

		% compute edges structure
		[obj.sigma,obj.edges,obj.mid_edges] = obj.str_sigma_2d(obj.nodes, obj.cells, obj.cc);
		obj.nsig = size(obj.sigma,1);
        obj.nedges = size(obj.edges,1);

		% compute indices 
		obj.ind = obj.indices(obj.ncells,obj.sigma);
		obj.nsig_in = obj.ind.nsig_in;

    end

		
		function [refined, coarse2fine] = refine(obj)

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

			% shorthands
			ncells = obj.ncells;
			nsig = obj.nsig;
			nnode = obj.nnodes;
			

			
			% compute cells-edges structure
			maxn= max(obj.cells(:,1));
			[cell_sig,~,~] = obj.str_cell_2d(maxn,obj.ind,obj.sigma,obj.mid_edges,obj.cc);
		
			cells_f = zeros(3*ncells,5);
			nodes_f = zeros(nnode+nsig+ncells,2);
			% TODO: remove this, derived variables are computed in the constructor
			nodes = [obj.nodes zeros(nnode,1)];
			mid = [obj.mid_edges zeros(nsig,1)];
			cc = [obj.cc zeros(ncells,1)];
			int = 0;
			coarse2fine = sparse(3*ncells,ncells);

			for k=1:ncells
				
				K = obj.cells(k,2:4);
				E = cell_sig(k,2:4);
				
				for j=1:3
					
					
					ej = E(obj.edges(E',1)==K(j) | obj.edges(E',2)==K(j));
					
					if nodes(K(j),3)==0
            int = int+1;
            A = int;
            nodes(K(j),3) = int;
            nodes_f(int,:) = nodes(K(j),1:2);
					else
            A = nodes(K(j),3);
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
			refined.init(cells_f(:,2:end),nodes_f);
			
		end
			
   	% destructor
		function obj = kill(obj)
			if (obj.initalized == 1)
				clear obj.cells;
				clear obj.edges;
                clear obj.sigma;
				clear obj.ind;
				clear obj.area;
				clear obj.cc;
				clear obj.mid;
			end
		end


		function [mass,edge_mass,div,grad] = build_matrices(obj)
            % create the matrices used in the finite volume discretization
		    % INPUTS
		    % obj : mesh
		    % OUTPUTS:
		    % mass : diagonal matrix with cell volume
		    % edge_mass : diagonal matrix with edge length
		    % div : sparse matrix discretizing the divergence operator
		    % grad : sparse matrix discretizing the grad operator
		    %        (grad = - edge_mass^{-1} div )
			mass = spdiags(obj.area,0,obj.ncells,obj.ncells);
			nsig_in = obj.nsig_in;
			ds = obj.sigma(obj.ind.internal,3).*obj.sigma(obj.ind.internal,4);
			edge_mass = spdiags(ds,0,nsig_in,nsig_in);
			div = Div2D(obj.ncells,nsig_in,obj.ind,obj.sigma); % divergence matrix
			grad = -edge_mass\div'; % gradient matrix
		end
		
		% info
		function obj = info(obj,fid)
            if (~exist('fid','var') )
				fid=1;
			end
		end

        function [cc,area,h] = mesh_2d(obj,nodes,cells)

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
            
            
            % % compute the orthogonal vector to the first edge vertex_1-vertex_2
            % v = [nodes(vert(:,1),2)-nodes(vert(:,2),2) nodes(vert(:,2),1)-nodes(vert(:,1),1)];
            % % use the second edge vertex_2-vertex_3 if a node's x-ccordinate is equal to zero
            % corrv = find(v(:,1)==0);
            % v(corrv,:) = [nodes(vert(corrv,2),2)-nodes(vert(corrv,3),2) nodes(vert(corrv,3),1)-nodes(vert(corrv,2),1)];
            % if any(v(:,1)==0)
            %     error('there is a problem with the mesh')
            % end
            % 
            % % compute the orthogonal vector to the last edge vertex_1-vertex_nvert
            % w = zeros(size(vert,1),2);
            % for i=1:size(vert,1)
            %     w(i,:) = [nodes(vert(i,1),2)-nodes(vert(i,nvert(i)),2) nodes(vert(i,nvert(i)),1)-nodes(vert(i,1),1)];
            % end
            % % detect parallel edges (v//w) and change the edge w: use instead the edge
            % % vertex_(nvert-1)-vertex_nvert 
            % corrw = [];
            % for i=1:size(v)
            %     if w(i,:)-v(i,:)==[0 0]
            %         corrw = [corrw; i];
            %     end
            % end
            % for i=1:length(corrw)
            %     w(corrw(i),:) = [nodes(vert(corrw(i),nvert(corrw(i))-1),2)-nodes(vert(corrw(i),nvert(corrw(i))),2) nodes(vert(corrw(i),nvert(corrw(i))),1)-nodes(vert(corrw(i),nvert(corrw(i))-1),1)];
            % end
            % 
            % % compute the midpoints of the two edges considered
            % midv = [0.5*(nodes(vert(:,1),1)+nodes(vert(:,2),1)) 0.5*(nodes(vert(:,2),2)+nodes(vert(:,1),2))];
            % midv(corrv,:) = [0.5*(nodes(vert(corrv,2),1)+nodes(vert(corrv,3),1)) 0.5*(nodes(vert(corrv,3),2)+nodes(vert(corrv,2),2))];
            % midw = zeros(size(vert,1),2);
            % for i=1:size(vert,1)
            %     midw(i,:) = [0.5*(nodes(vert(i,1),1)+nodes(vert(i,nvert(i)),1)) 0.5*(nodes(vert(i,nvert(i)),2)+nodes(vert(i,1),2))];
            % end
            % for i=1:length(corrw)
            %     midw(corrw(i),:) = [0.5*(nodes(vert(corrw(i),nvert(corrw(i))-1),1)+nodes(vert(corrw(i),nvert(corrw(i))),1)) 0.5*(nodes(vert(corrw(i),nvert(corrw(i))),2)+nodes(vert(corrw(i),nvert(corrw(i))-1),2))];
            % end
            % 
            % % compute the point of intersection between the two vectors
            % beta = (midw(:,2)-midv(:,2)-v(:,2).*(midw(:,1)-midv(:,1))./v(:,1))./(v(:,2).*w(:,1)./v(:,1)-w(:,2));
            % cc = midw + beta.*w;
            
            % formula for the circumcenter of a triangle
            % (for general regular polygons we take the first three nodes to find the
            % circumscribed circle)
            a = [nodes(cells(:,2),1) nodes(cells(:,2),2) zeros(size(cells,1),1)];
            b = [nodes(cells(:,3),1) nodes(cells(:,3),2) zeros(size(cells,1),1)];
            c = [nodes(cells(:,4),1) nodes(cells(:,4),2) zeros(size(cells,1),1)];
            % a = [nodes(cells(:,1),1) nodes(cells(:,1),2) nodes(cells(:,1),3)];
            % b = [nodes(cells(:,2),1) nodes(cells(:,2),2) nodes(cells(:,2),3)];
            % c = [nodes(cells(:,3),1) nodes(cells(:,3),2) nodes(cells(:,3),3)];
            ac = c-a;
            ab = b-a;
            abxac = cross(ab,ac,2);
            abxacxab = cross(abxac,ab,2);
            acxabxac = cross(ac,abxac,2);
            cc = zeros(size(cells,1),3);
            for k=1:size(cells,1)
                cc(k,:) = a(k,:) + (norm(ac(k,:),2)^2*abxacxab(k,:) + norm(ab(k,:),2)^2*acxabxac(k,:))/(2*norm(abxac(k,:),2)^2);
            end
            cc=cc(:,1:2);
            
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

        function [sigma,edges,mid_edges] = str_sigma_2d(obj,nodes,cells,cc)

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
            mid_edges = [0.5*(nodes(edges(:,1),1)+nodes(edges(:,2),1)) 0.5*(nodes(edges(:,1),2)+nodes(edges(:,2),2))];
            
            
            % add the edges' distances d_sig, length m_sig and the transmissivity d_sig/m_sig
            % add also the distances d_{K,sig}, d_{L,sig}
            edges = [edges zeros(size(edges,1),5)];
            for sig=1:size(edges,1)
                if edges(sig,4)~=0
                    d_sig = sqrt( (cc(edges(sig,3),1)-cc(edges(sig,4),1)).^2 + (cc(edges(sig,3),2)-cc(edges(sig,4),2)).^2 );
                    dKsig = sqrt((mid_edges(sig,1)-cc(edges(sig,3),1)).^2+...
                        (mid_edges(sig,2)-cc(edges(sig,3),2)).^2);
                    dLsig = sqrt((mid_edges(sig,1)-cc(edges(sig,4),1)).^2+...
                        (mid_edges(sig,2)-cc(edges(sig,4),2)).^2);
                else
                    d_sig = sqrt( (cc(edges(sig,3),1)-mid_edges(sig,1)).^2 + (cc(edges(sig,3),2)-mid_edges(sig,2)).^2 );
                    dKsig = sqrt((mid_edges(sig,1)-cc(edges(sig,3),1)).^2+...
                        (mid_edges(sig,2)-cc(edges(sig,3),2)).^2);
                    dLsig = 0;
                end
                m_sig = sqrt((nodes(edges(sig,1),1)-nodes(edges(sig,2),1)).^2+(nodes(edges(sig,1),2)-nodes(edges(sig,2),2)).^2);
                edges(sig,5:9) = [d_sig m_sig dKsig dLsig m_sig/d_sig];
            end
            
            sigma = edges(:,3:end);
            edges = edges(:,1:2);

		end

		function [ind] = indices(obj,ncells,sigma)
        
            % INPUT:
            % ncells: number of cells
            % edges: edges structure
            % OUPUT:
            % ind: indices structure
            
            % this function creates indices for a simplified access to the matrices
            % these indices enable to select the entries of submatrices of the FV
            % matrices
            
            % indices for internal edges:
            
            internal = find(sigma(:,2)~=0);
            
            ind.internal = internal;
            ind.bound = find(sigma(:,2)==0);
            
            ind.nsig_in = length(internal);
            ind.nsig_b = length(ind.bound);
            
            % OLD
            
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

		
		function [cell_sig,cell_sig_in,cell_dist] = str_cell_2d(obj,maxn,ind,sigma,mid,cc)

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
            
            ncells = size(cc,1);
            nsigma = size(sigma,1);
            nsig_in = length(ind.internal);
            
            cell_sig_in = zeros(ncells,maxn+1);
            % for each internal edge e, take the two cells K and L defining it and assign to both the
            % edge e
            for sig=1:nsig_in
                K = sigma(ind.internal(sig),1);
                L = sigma(ind.internal(sig),2);
                cell_sig_in(K,1) = cell_sig_in(K,1)+1;
                cell_sig_in(K,cell_sig_in(K,1)+1) = sig;
                if L~=0
                    cell_sig_in(L,1) = cell_sig_in(L,1)+1;
                    cell_sig_in(L,cell_sig_in(L,1)+1) = sig;
                end
            end
            
            cell_sig = zeros(ncells,maxn+1);
            % for each edge sig, take the two cells K and L defining it and assign to both the
            % edge sig
            for sig=1:nsigma
                K = sigma(sig,1);
                L = sigma(sig,2);
                cell_sig(K,1) = cell_sig(K,1)+1;
                cell_sig(K,cell_sig(K,1)+1) = sig;
                if L~=0
                    cell_sig(L,1) = cell_sig(L,1)+1;
                    cell_sig(L,cell_sig(L,1)+1) = sig;
                end
            end
            
            cell_dist = zeros(ncells,maxn+1);
            cell_dist(:,1) = cell_sig(:,1);
            % for each cell i, for each edges j of the cell i, compute the distance |cc(i)-mid(j)|
            for i=1:ncells
                for j=1:cell_dist(i,1)
                    cell_dist(i,1+j) = sqrt((cc(i,1)-mid(cell_sig(i,1+j),1)).^2+(cc(i,2)-mid(cell_sig(i,1+j),2)).^2);
                end
            end

        end
		
	end
end

function [area] = area_pol(Np,coord)

%compute the area of the polygon
vertices = [coord; coord(1,:)];
sum1=0; sum2=0;
for i=1:Np
    sum1 = sum1 + vertices(i,1)*vertices(i+1,2);
    sum2 = sum2 + vertices(i,2)*vertices(i+1,1);
end
area = (sum1-sum2)/2;
area = abs(area);

end
