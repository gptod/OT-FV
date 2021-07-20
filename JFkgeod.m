function [JOC] = JFkgeod(N,Dt,divt,Mxt,Mxt2h,gradt,It,rhos,drhos,ddrhosa,uk)

% Compute jacobian of the system of first-order optimality conditions

Nt = N+1;
ncell2h = size(Mxt2h,1)/Nt;
ncell = size(Mxt,1)/Nt;
nei = size(divt,2)/Nt;
tnr2h = N*ncell2h;
tnp = Nt*ncell;
te = Nt*nei;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);
s = uk(tnp+tnr2h+1:tnp+2*tnr2h);

I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
gradphi = gradt*phi;

% derivatives of the optimality conditions in phi
JOCpp = divt*spdiags(rhos,0,te,te)*gradt;
JOCpr = Mxt*Nt*Dt*It*I_all+divt*spdiags(gradphi,0,te,te)*drhos;
JOCps = sparse(tnp,tnr2h);

% derivatives of the optimality conditions in rho
%JOCrp = (Nt*Mxt*Dt*It*I_all)'-drhos'*Mst*(spdiags(gradphi,0,te,te)*gradt);
JOCrr = -0.5*ddrhosa;
JOCrs = -Mxt2h(1:tnr2h,1:tnr2h);

% optimality conditions in s
JOCsp = sparse(tnr2h,tnp);
JOCsr = -spdiags(s,0,tnr2h,tnr2h);
JOCss = -spdiags(rho,0,tnr2h,tnr2h);



JOC.pp = JOCpp; JOC.pr = JOCpr; JOC.ps = JOCps;
JOC.rp = JOCpr'; JOC.rr = JOCrr; JOC.rs = JOCrs;
JOC.sp = JOCsp; JOC.sr = JOCsr; JOC.ss = JOCss;

JOC.ntimestep=N;
JOC.ncellrho = ncell2h;
JOC.ncellphi = ncell;

JOC.DtT=Nt*Dt*It*I_all;
JOC.B1T_time=Mxt*Nt*Dt*It*I_all;
JOC.B1T_space=divt*spdiags(gradphi,0,te,te)*drhos;
  
JOC.LT=spdiags(gradphi,0,te,te)*drhos;
JOC.diagrho=spdiags(rhos,0,te,te);
% for i=2:N-1
%   irow=i
%   jcol=i-1
%   diff=JOC.LT(1+(irow-1)*ncell:(irow+1)*ncell,1+(jcol-1)*ncell2h:(jcol)*ncell2h)-...
%        JOC.LT(1+(irow-1)*ncell:(irow+1)*ncell,1+(jcol-1)*ncell2h:(jcol)*ncell2h);
%   norm(full(diff))
% end

JOC.It=It;

JOC.Mxt=Mxt;


end
