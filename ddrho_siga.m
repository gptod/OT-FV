function [drhosa] = ddrho_siga(ind,edges,mid,cc,rho,a)

% second derivative of the harmonic reconstruction from cells to edges applied to
% the vector a

ncell = size(cc,1);
nei = length(ind.internal);
drhosa = sparse(ncell,ncell);

for e=1:nei
    
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    dK = sqrt((mid(ind.internal(e),1)-cc(K,1))^2+(mid(ind.internal(e),2)-cc(K,2))^2);
    dL = sqrt((mid(ind.internal(e),1)-cc(L,1))^2+(mid(ind.internal(e),2)-cc(L,2))^2);
    ds = edges(ind.internal(e),5);
    % harmonic mean
    drhosa(K,K) = drhosa(K,K) + dKdKharm(rho(K),rho(L),dK,dL,ds)*a(e);
    drhosa(K,L) = drhosa(K,L) + dLdKharm(rho(K),rho(L),dK,dL,ds)*a(e);
    drhosa(L,K) = drhosa(L,K) + dLdKharm(rho(K),rho(L),dK,dL,ds)*a(e);
    drhosa(L,L) = drhosa(L,L) + dKdKharm(rho(L),rho(K),dL,dK,ds)*a(e);

end




