function RHt = assembleRHt(N,ncell)

% time reconstruction matrix 

RHt = sparse((N+1)*ncell,(N+2)*ncell);

for k=1:N+1
   
    % arithmetic mean as reconstruction at the staggered time 
    RHt((k-1)*ncell+1:k*ncell,(k-1)*ncell+1:(k+1)*ncell) = [0.5*speye(ncell,ncell) 0.5*speye(ncell,ncell)];
    % final/initial density as reconstruction at the staggered time
    %RHt((k-1)*ncell+1:k*ncell,(k-1)*ncell+1:(k+1)*ncell) = [sparse(ncell,ncell) speye(ncell,ncell)];
    %RHt((k-1)*ncell+1:k*ncell,(k-1)*ncell+1:(k+1)*ncell) = [speye(ncell,ncekk) sparse(ncell,ncell)];

end
