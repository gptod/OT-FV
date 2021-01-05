function Mst = assembleMst(N,nedge,Ms)


Mst = sparse((N+1)*nedge,(N+1)*nedge);

for k=1:N+1
   
    Mst((k-1)*nedge+1:k*nedge,(k-1)*nedge+1:k*nedge) = Ms;
    
end