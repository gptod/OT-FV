function [dnodes,dcells] = diamond_2d(nodes,cc,edges,sigma)

nnode = size(nodes,1);
nsig = size(sigma,1);
dnodes = [nodes; cc];
dcells = zeros(nedge,5);

for sig=1:nsig
    
    if sigma(sig,2)==0
        A = edges(sig,1);
        B = edges(sig,2);
        C = nnode+sigma(sig,1);
        D = A;
        dcells(sig,:) = [3 A C B D];
    else
        A = edges(sig,1);
        B = edges(sig,2);
        C = nnode+sigma(sig,1);
        D = nnode+sigma(sig,2);
        dcells(sig,:) = [4 A C B D];
    end
    
end