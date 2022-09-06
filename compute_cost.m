function cost = compute_cost(gradt,Mst,N,rhos,phi)

% Compute the Wasserstein distance, namely the squared root of twice the
% approximated transport cost

Nt = N+1;
gradphi = gradt*phi;
cost = (1/Nt)*rhos'*Mst*(gradphi.^2); cost = sqrt(cost);