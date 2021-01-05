function Dt = assembleDt(N,ncell)


Dt = sparse((N+1)*ncell,(N+2)*ncell);


for k=1:N+1
   
    Dt((k-1)*ncell+1:k*ncell,(k-1)*ncell+1:(k+1)*ncell) = [-speye(ncell,ncell) speye(ncell,ncell)];
    
end