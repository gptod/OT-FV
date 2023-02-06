function A = FVstiff(ind,sigma,ncell)

KK=sigma(ind.internal,1);
LL=sigma(ind.internal,2);

HKK = sigma(ind.internal,7);
HLL = sigma(ind.internal,7);
HKL = -sigma(ind.internal,7);
HLK = -sigma(ind.internal,7);


A = sparse(KK,KK,HKK,ncell,ncell);
A = A + sparse(KK,LL,HKL,ncell,ncell);
A = A + sparse(LL,KK,HLK,ncell,ncell);
A = A + sparse(LL,LL,HLL,ncell,ncell);
