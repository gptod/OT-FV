clear
close all


load('tri2_mesh1.mat'); % Delaunay mesh
clear name fn1 fn2

name = 'vor_mesh1'; % name of the final structure

% specify the corners of the mesh:
bl = find(nodes(:,1)==min(nodes(:,1))&nodes(:,2)==min(nodes(:,2)));
br = find(nodes(:,1)==max(nodes(:,1))&nodes(:,2)==min(nodes(:,2)));
tr = find(nodes(:,1)==max(nodes(:,1))&nodes(:,2)==max(nodes(:,2)));
tl = find(nodes(:,1)==min(nodes(:,1))&nodes(:,2)==max(nodes(:,2)));
corners = [bl; br; tr; tl];

% shift corners and boundary midpoints to avoid singular diamond cells
sh = 0.2*h; %length of the shift (h meshsize)
for e=1:length(ind.bound)
    ee = ind.bound(e);
    % compute the orthogonal vector to the edge
    v = [nodes(edges(ee,1),2)-nodes(edges(ee,2),2) nodes(edges(ee,2),1)-nodes(edges(ee,1),1)];
    v = -v*sum(v.*(cc(edges(ee,3),:)-mid(ee,:)));
    v = v/norm(v);
    % shift the midpoint of l towards the outside 
    mid(ee,:) = mid(ee,:)+sh*v;
    if any(nodes(edges(ee,1),:)==nodes(corners,:))
        nodes(edges(ee,1),:) = nodes(edges(ee,1),:)+sh*v;
    elseif any(nodes(edges(ee,2),:)==nodes(corners,:))
        nodes(edges(ee,2),:) = nodes(edges(ee,2),:)+sh*v;
    end
end

% generate a voronoi diagram out of the Delaunay triangulation
[nodes_v,cells_v] = voronoi(nodes,edges,cc,mid,corners);
% rescale all the coordinates
nodes_v(:,1) = (nodes_v(:,1)+sh)/(1+2*sh);
nodes_v(:,2) = (nodes_v(:,2)+sh)/(1+2*sh);
% clear all the old variables
clear e ee v sh corners bl br tr tl
clear cell_dist cell_e cell_eint cells mid cc dcells dnodes edges
clear h ind maxn ncell nedge nei nnode nodes area fn1

nodes = nodes_v; cells = cells_v;
clear nodes_v cells_v
nnode = size(nodes,1);
ncell = size(cells,1);
maxn= max(cells(:,1));

% compute cells related properties
[cc,area,h] = mesh(nodes,cells);

% compute edges structure
[edges,mid] = str_edge(nodes,cells,cc);
nedge = size(edges,1);

% compute indices
ind = indices(ncell,edges);

% compute cells-edges structure
nei = length(ind.internal);
[cell_e,cell_eint,cell_dist] = str_cell(maxn,ind,edges,mid,cc);

% compute diamond cells structure
[dnodes,dcells] = diamond(nnode,nedge,nodes,cc,edges);

% print mesh
patch('Faces',edges(:,1:2),'Vertices',nodes)
axis square
axis off

% save the mesh structure
save(name)



