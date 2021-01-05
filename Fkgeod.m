function OC = Fkgeod(ind,edges,cc,mid,N,rho_f,rho_in,Dt,divt,Mxt,Mxt2h,Mst,gradt,RHt,It,rec,uk,mu)

% Compute the system of first-order optimality conditions

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
else
    rhoa = RHt*It*rho_all;
    rhos = zeros(te,1);
    drhos = sparse(te,tnr);
    for k=1:N+1
        rhos((k-1)*nei+1:k*nei) = rho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
        drhos((k-1)*nei+1:k*nei,(k-1)*ncell+1:k*ncell) = drho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
    end
    drhos = drhos*RHt*It*I_all;
end


% optimality conditions in phi
OCp = Mxt*Nt*Dt*It*rho_all + divt*(rhos.*gradphi);

% optimality conditions in rho
OCr = (Nt*Mxt*Dt*It*I_all)'*phi - (0.5*drhos'*Mst*(gradphi).^2) - Mxt2h(1:tnr2h,1:tnr2h)*s;

% optimality conditions in s
OCs = (-rho.*s+mu);


OC.p = OCp;
OC.r = OCr;
OC.s = OCs;


end

