function [cell_sig,cell_sig_in,cell_dist] = str_cell_2d(maxn,ind,sigma,mid,cc)

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
nsigma = size(sigma,1);
nsig_in = length(ind.internal);

cell_sig_in = zeros(ncell,maxn+1);
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

cell_sig = zeros(ncell,maxn+1);
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


cell_dist = zeros(ncell,maxn+1);
cell_dist(:,1) = cell_sig(:,1);
% for each cell i, for each edges j of the cell i, compute the distance |cc(i)-mid(j)|
for i=1:ncell
    for j=1:cell_dist(i,1)
        cell_dist(i,1+j) = sqrt((cc(i,1)-mid(cell_sig(i,1+j),1)).^2+(cc(i,2)-mid(cell_sig(i,1+j),2)).^2);
    end
end
