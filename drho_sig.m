function [drhos] = drho_sig(ind,sigma,rho)

% derivative of the harmonic reconstruction from cells to edges


ncell = size(rho,1);
nsig_in = length(ind.internal);

KK=sigma(ind.internal,3);
LL=sigma(ind.internal,4);

dK = sigma(ind.internal,5);
dL = sigma(ind.internal,6);
ds = sigma(ind.internal,3);


% drhos2 = sparse(nsig_in,ncell);
% drhos2(ind.et_K) = dKharm(rho(KK),rho(LL),dK,dL,ds);
% drhos2(ind.et_L) = dKharm(rho(LL),rho(KK),dL,dK,ds);

% this assigning is faster than the other one
drhos = sparse(1:nsig_in,KK,dKharm(rho(KK),rho(LL),dK,dL,ds),nsig_in,ncell);
drhos = drhos + sparse(1:nsig_in,LL,dKharm(rho(LL),rho(KK),dL,dK,ds),nsig_in,ncell);



