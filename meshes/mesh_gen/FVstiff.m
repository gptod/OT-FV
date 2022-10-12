function A = FVstiff(ind,edges,ncell)

KK=edges(ind.internal,3);
LL=edges(ind.internal,4);

HKK = edges(ind.internal,7);
HLL = edges(ind.internal,7);
HKL = -edges(ind.internal,7);
HLK = -edges(ind.internal,7);


A = sparse(KK,KK,HKK,ncell,ncell);
A = A + sparse(KK,LL,HKL,ncell,ncell);
A = A + sparse(LL,KK,HLK,ncell,ncell);
A = A + sparse(LL,LL,HLL,ncell,ncell);
