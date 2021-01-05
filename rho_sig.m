function [rhos] = rho_sig(ind,edges,mid,cc,rho)

% harmonic reconstruction from cells to edges

nei = length(ind.internal);
rhos = zeros(nei,1);

for e=1:nei
    
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    dK = sqrt((mid(ind.internal(e),1)-cc(K,1))^2+(mid(ind.internal(e),2)-cc(K,2))^2);
    dL = sqrt((mid(ind.internal(e),1)-cc(L,1))^2+(mid(ind.internal(e),2)-cc(L,2))^2);
    ds = edges(ind.internal(e),5);
    % harmonic mean
    rhos(e) = rhos(e) + harm(rho(K),rho(L),dK,dL,ds);

end


