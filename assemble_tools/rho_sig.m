function [rhos] = rho_sig(ind,sigma,rho)

% harmonic reconstruction from cells to edges

KK=sigma(ind.internal,1);
LL=sigma(ind.internal,2);

dK = sigma(ind.internal,5);
dL = sigma(ind.internal,6);
ds = sigma(ind.internal,3);

rhos = harm(rho(KK),rho(LL),dK,dL,ds);
