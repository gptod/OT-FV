function [dnodes,dcells] = diamond(nnode,nedge,nodes,cc,edges,mid)

dnodes = [nodes; cc];
dcells = zeros(nedge,5);

for e=1:nedge
    
    if edges(e,4)==0
        A = edges(e,1);
        B = edges(e,2);
        C = nnode+edges(e,3);
        D = A;
        dcells(e,:) = [3 A C B D];
    else
        A = edges(e,1);
        B = edges(e,2);
        C = nnode+edges(e,3);
        D = nnode+edges(e,4);
        dcells(e,:) = [4 A C B D];
    end
    
end