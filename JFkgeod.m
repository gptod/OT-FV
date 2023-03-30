function [JOC] = JFkgeod(N,Dt,divt,Mxt,Mxt2h,gradt,It,rhos,drhos,ddrhosa,uk,I,Rs,Rst_rho,divt_rho,gradt_rho)

% Compute jacobian of the system of first-order optimality conditions

Nt = N+1;
ncell2h = size(Mxt2h,1)/Nt;
ncell = size(Mxt,1)/Nt;
nsig_in = size(divt,2)/Nt;
tnr2h = N*ncell2h;
tnp = Nt*ncell;
te = Nt*nsig_in;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);
s = uk(tnp+tnr2h+1:tnp+2*tnr2h);


I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
gradphi = gradt*phi;

% derivatives of the optimality conditions in phi
JOCpp = divt*spdiags(rhos,0,te,te)*gradt;
JOCpr = Mxt*Nt*Dt*It*I_all+divt*spdiags(gradphi,0,te,te)*drhos; %drhos=Rst*RHt*It*I_all
JOC.Dt = Mxt*Nt*Dt*It*I_all;
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
JOC.G=spdiags(gradphi,0,te,te)*drhos;

%temp=(spdiags(gradphi,0,te,te)*drhos)'*spdiags(gradphi,0,te,te)*drhos;
%imagesc(log(temp))
%return


JOC.LT=spdiags(gradphi,0,te,te)*drhos;
JOC.diagrho=spdiags(rhos,0,te,te);


JOC.It=It;

JOC.Mxt=Mxt;


				% store this diagonal
Np=tnp;
JOC.area2h=-spdiags(JOC.rs(1:ncell2h,1:ncell2h));
JOC.area=spdiags(Mxt(1:ncell,1:ncell));
JOC.areadomain=sum(JOC.area2h);



JOC.I=I;

JOC.div=divt(1:ncell,1:size(divt,2)/Nt);
JOC.grad=gradt(1:size(divt,2)/Nt,1:ncell);
JOC.Rst_rho=Rst_rho;
JOC.Rs=Rs;
JOC.divt_rho=divt_rho;
JOC.gradt_rho=gradt_rho;

JOC.divt_phi=divt;
JOC.gradt_phi=gradt;



JOC.gradphi=gradphi;
JOC.phi = phi;
end
