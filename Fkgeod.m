function OC = Fkgeod(N,rho_f,rho_in,Dt,divt,Mxt,Mxt2h,Mst,gradt,It,rhosk,drhosk,uk,mu)
  
% Compute the system of first-order optimality conditions
Nt = N+1;
ncell2h = size(Mxt2h,1)/Nt;
ncell = size(Mxt,1)/Nt;
tnr2h = N*ncell2h;
tnp = Nt*ncell;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);
s = uk(tnp+tnr2h+1:tnp+2*tnr2h);

rho_all = [rho_in;rho;rho_f];
I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
gradphi = gradt*phi;

% optimality conditions in phi
OCp = Mxt*Nt*Dt*It*rho_all + divt*(rhosk.*gradphi);

% optimality conditions in rho
OCr = (Nt*Mxt*Dt*It*I_all)'*phi - (0.5*drhosk'*Mst*(gradphi).^2) - Mxt2h(1:tnr2h,1:tnr2h)*s;

% optimality conditions in s
OCs = (-rho.*s+mu);

% create data structure
OC.p = OCp;
OC.r = OCr;
OC.s = OCs;

% copy sizes
OC.N = N;
OC.ncell2h=ncell2h;
OC.ncell=ncell;

end

