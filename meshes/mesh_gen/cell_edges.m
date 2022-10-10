function [cell_e] = cell_edges(ntri,nedge,edges)


cell_e = zeros(ntri,4);
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
