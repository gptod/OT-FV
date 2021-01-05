function Mxt = assembleMxt(N,ncell,Mx)


Mxt = sparse((N+1)*ncell,(N+1)*ncell);

for k=1:N+1
   
    Mxt((k-1)*ncell+1:k*ncell,(k-1)*ncell+1:k*ncell) = Mx;
    
end