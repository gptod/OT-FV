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
