function [drhos] = drho_sig(ind,edges,mid,cc,rho)

% derivative of the harmonic reconstruction from cells to edges


ncell = size(cc,1);
nei = length(ind.internal);
drhos = sparse(nei,ncell);

for e=1:nei
    
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    dK = sqrt((mid(ind.internal(e),1)-cc(K,1))^2+(mid(ind.internal(e),2)-cc(K,2))^2);
    dL = sqrt((mid(ind.internal(e),1)-cc(L,1))^2+(mid(ind.internal(e),2)-cc(L,2))^2);
    ds = edges(ind.internal(e),5);
    % harmonic mean
    drhos(e,K) = drhos(e,K) + dKharm(rho(K),rho(L),dK,dL,ds);
    drhos(e,L) = drhos(e,L) + dKharm(rho(L),rho(K),dL,dK,ds);
        
end


