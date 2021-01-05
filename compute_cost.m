function cost = compute_cost(ind,edges,mid,cc,gradt,Mst,RHt,St,N,rho_all,phi,rec)

% Compute the Wasserstein distance, namely the squared root of twice the
% approximated transport cost

ncell = size(cc,1);
nei = length(ind.internal);
Nt = N+1;

gradphi = gradt*phi;
if rec==1
    Rst = sparse((N+1)*nei,(N+1)*ncell);
    for k=1:N+1
        Rst((k-1)*nei+1:k*nei,(k-1)*ncell+1:k*ncell) = Ktos2D(ind,edges,cc,mid);
    end
    rhoa = RHt*St*rho_all;
    rhos = Rst*rhoa;
else
    rhoa = RHt*St*rho_all;
    rhos = zeros((N+1)*nei,1);
    for k=1:N+1
        rhos((k-1)*nei+1:k*nei) = rho_sig(ind,edges,mid,cc,rhoa((k-1)*ncell+1:k*ncell));
    end
end

cost = (1/Nt)*rhos'*Mst*(gradphi.^2); cost = sqrt(cost);