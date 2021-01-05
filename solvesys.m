function d = solvesys(JF,F,indc)

% solve the linear system d = -JF\F
% impose the condition d(indc) = 0 to fix the constant of the potential dp

Np = length(F.p);
Nr = length(F.r);
Ns = length(F.s);

S = [JF.pp JF.pr; JF.rp JF.rr-JF.rs*(JF.ss\JF.sr)];
S(indc,:) = sparse(1,Np+Nr); S(:,indc) = sparse(Np+Nr,1); S(indc,indc) = 1;
f1 = F.p; f1(indc) = 0;
f2 = F.r-JF.rs*(JF.ss\F.s);
d = S\[-f1;-f2];
dp = d(1:Np); dr = d(Np+1:end);
ds = JF.ss\(-F.s-JF.sr*dr);
d = [dp; dr; ds];


