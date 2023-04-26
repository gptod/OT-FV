function [Rs] = Ktos2D(ind,sigma,cc,mid)

% reconstruction operator from cells to edges

ncells = size(cc,1);
nsig_in = length(ind.internal);
Rs = sparse(nsig_in,ncells);

for sig=1:nsig_in
    K = sigma(ind.internal(sig),1);
    L = sigma(ind.internal(sig),2);
    % mass average
    %dK = sqrt((mid(ind.internal(sig),1)-cc(K,1))^2+(mid(ind.internal(sig),2)-cc(K,2))^2);
    %dL = sqrt((mid(ind.internal(sig),1)-cc(L,1))^2+(mid(ind.internal(sig),2)-cc(L,2))^2);
    dK = sigma(ind.internal(sig),5);
    dL = sigma(ind.internal(sig),6);
    ds = sigma(ind.internal(sig),3);
    Rs(sig,K) = Rs(sig,K) + dK/ds;
    Rs(sig,L) = Rs(sig,L) + dL/ds;
    %Rs(e,K) = Rs(e,K) + dL/ds;
    %Rs(e,L) = Rs(e,L) + dK/ds;
    %Rs(e,K) = Rs(e,K) + 0.5;
    %Rs(e,L) = Rs(e,L) + 0.5;
end


