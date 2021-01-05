function [Rs] = Ktos2D(ind,edges,cc,mid)

% reconstruction operator from cells to edges

ncell = size(cc,1);
nei = length(ind.internal);
Rs = sparse(nei,ncell);

for e=1:nei
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    % mass average
    dK = sqrt((mid(ind.internal(e),1)-cc(K,1))^2+(mid(ind.internal(e),2)-cc(K,2))^2);
    dL = sqrt((mid(ind.internal(e),1)-cc(L,1))^2+(mid(ind.internal(e),2)-cc(L,2))^2);
    ds = edges(ind.internal(e),5);
    Rs(e,K) = Rs(e,K) + dK/ds;
    Rs(e,L) = Rs(e,L) + dL/ds;
    %Rs(e,K) = Rs(e,K) + dL/ds;
    %Rs(e,L) = Rs(e,L) + dK/ds;
    %Rs(e,K) = Rs(e,K) + 0.5;
    %Rs(e,L) = Rs(e,L) + 0.5;
end


