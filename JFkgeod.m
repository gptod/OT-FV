function [JOC] = JFkgeod(ind,edges,cc,mid,N,rho_f,rho_in,Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk)

% Compute jacobian of the system of first-order optimality conditions

Nt = N+1;
ncell2h = size(rho_in,1);
ncell = size(cc,1);
nei = length(ind.internal);
tnr = N*ncell;
tnr2h = N*ncell2h;
tnp = (N+1)*ncell;
te = (N+1)*nei;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);
s = uk(tnp+tnr2h+1:tnp+2*tnr2h);

rho_all = [rho_in;rho;rho_f];
I_all = [zeros(ncell2h,tnr2h);speye(tnr2h,tnr2h);zeros(ncell2h,tnr2h)];
gradphi = gradt*phi;

% compute the reconstructed density on the diamond cells, rhos
if rec==1
    Rst = sparse(te,tnp);
    for k=1:N+1
        Rst((k-1)*nei+1:k*nei,(k-1)*ncell+1:k*ncell) = Ktos2D(ind,edges,cc,mid);
    end    
    rhos = Rst*RHt*It*rho_all;
    drhos = Rst*RHt*It*I_all;
    ddrhosa = sparse(tnr2h,tnr2h);
else
    rhoa = RHt*It*rho_all;
    rhos = zeros(te,1);
    drhos = sparse(te,tnr);
    ddrhosa = sparse(tnr,tnr);
    a = Mst*gradphi.^2;
    for k=1:N+1
        rhos((k-1)*nei+1:k*nei) = rho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
        drhos((k-1)*nei+1:k*nei,(k-1)*ncell+1:k*ncell) = ...
            drho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
        ddrhosa((k-1)*ncell+1:k*ncell,(k-1)*ncell+1:k*ncell) = ...
            ddrho_siga(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell),a((k-1)*nei+1:k*nei));
    end
    drhos = drhos*RHt*It*I_all;
    ddrhosa = I_all'*It'*RHt'*ddrhosa*RHt*It*I_all;
end


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
JOC.ncellrho = size(rho_in,1);
JOC.ncellphi = size(cc,1);

end
