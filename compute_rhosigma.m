function [out]=compute_rhosigma(ind,sigma,cc,N,rho_f,rho_in,gradt,Mst,RHt,It,Rst,rec,uk,str)

% str='rhos'    -> compute reconstructed density
% str='drhos'   -> compute derivative of the reconstructed density
% str='ddrhosa' -> compute second derivative of the reconstructed density
%                  applied to the vector a

ncell2h = size(rho_in,1);
ncell = size(cc,1);
nsig_in = length(ind.internal);
tnr = N*ncell;
tnr2h = N*ncell2h;
tnp = (N+1)*ncell;
te = (N+1)*nsig_in;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);

rho_all = [rho_in;rho;rho_f];
I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
gradphi = gradt*phi;

% compute the reconstructed density on the diamond cells, rhos
if rec==1
    if strcmp(str,'rhos')
        out = Rst*RHt*It*rho_all;
    elseif strcmp(str,'drhos')
        out = Rst*RHt*It*I_all;
    elseif strcmp(str,'ddrhosa')
        out = sparse(tnr2h,tnr2h);
    end
else
    rhoa = RHt*It*rho_all;
    if strcmp(str,'rhos')
        out = zeros(te,1);
        for k=1:N+1
            out((k-1)*nsig_in+1:k*nsig_in) = rho_sig(ind,sigma,rhoa((k-1)*ncell+1:k*ncell));            
        end
    elseif strcmp(str,'drhos')
        out = cell(N+1,1);
        for k=1:N+1
            out{k} = drho_sig(ind,sigma,rhoa((k-1)*ncell+1:k*ncell));
        end
        out = blkdiag(out{:});
        out = out*RHt*It*I_all;
    elseif strcmp(str,'ddrhosa')
        out = cell(N+1,1);
        a = Mst*gradphi.^2;
        for k=1:N+1
            out{k} = ddrho_siga(ind,sigma,rhoa((k-1)*ncell+1:k*ncell),a((k-1)*nsig_in+1:k*nsig_in));
        end
        out = blkdiag(out{:});
        out = I_all'*It'*RHt'*out*RHt*It*I_all;
    end
end


