function [rhos,drhos,ddrhosa]=compute_rhosigma(ind,edges,cc,mid,N,rho_f,rho_in,gradt,Mst,RHt,It,Rst,rec,uk,str)

% Compute jacobian of the system of first-order optimality conditions

ncell2h = size(rho_in,1);
ncell = size(cc,1);
nei = length(ind.internal);
tnr = N*ncell;
tnr2h = N*ncell2h;
tnp = (N+1)*ncell;
te = (N+1)*nei;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);

rho_all = [rho_in;rho;rho_f];
I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
gradphi = gradt*phi;

% compute the reconstructed density on the diamond cells, rhos
if rec==1
    if strcmp(str,'rhos')
        rhos = Rst*RHt*It*rho_all;
        drhos = Rst*RHt*It*I_all;
        ddrhosa = [];
    elseif strcmp(str,'ddrhosa')
        rhos = [];
        drhos = [];
        ddrhosa = sparse(tnr2h,tnr2h);
    end
else
    rhoa = RHt*It*rho_all;
    if strcmp(str,'rhos')
        rhos = zeros(te,1);
        drhos = sparse(te,tnr);
        ddrhosa = [];
        for k=1:N+1
            rhos((k-1)*nei+1:k*nei) = rho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
            drhos((k-1)*nei+1:k*nei,(k-1)*ncell+1:k*ncell) = ...
                drho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
        end
        drhos = drhos*RHt*It*I_all;
    elseif strcmp(str,'ddrhosa')
        rhos = [];
        drhos = [];  
        ddrhosa = sparse(tnr,tnr);
        a = Mst*gradphi.^2;
        for k=1:N+1
            ddrhosa((k-1)*ncell+1:k*ncell,(k-1)*ncell+1:k*ncell) = ...
                ddrho_siga(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell),a((k-1)*nei+1:k*nei));
        end
        ddrhosa = I_all'*It'*RHt'*ddrhosa*RHt*It*I_all;
    end
end


