function [edges,mid] = str_edge(nodes,cells,cc)

% rough computation of the edges structure
% INPUT:
% node: nodes' coordinates
% cell: topology
% cc: cell orthocenters
% OUTPUT:
% edges: edges structure
% mid: midpoints vector

vert = cells(:,2:end); % disregard the column specifying the number of vertices per cell

nvert = cells(:,1);
nnode = size(nodes,1);
     
% construct all the edges with repetitions
% the edge is construct as: [vertexA vertexB], vertexA and vertexB are vertices numbers,
% with vertexA < vertexB
edges_rep = zeros(sum(nvert),3);
k = 0;
for i=1:size(vert,1)
    for j=1:nvert(i)
        k = k+1;
        if j==nvert(i)
            edges_rep(k,:) = [vert(i,j) vert(i,1) i];
        else
            edges_rep(k,:) = [vert(i,j) vert(i,j+1) i];
        end
    end
end
    
% order each edge from the smaller vertex to the bigger one
[edges_rep(:,1:2),~] = sort(edges_rep(:,1:2),2);
% order the edges in ascending order with respect to their first vertex
[~,ord2] = sort(edges_rep(:,1),1);
edges_rep = edges_rep(ord2,:);

% re-order the edges in ascending order also with respect to their second vertex
row = 1;
for i=1:nnode
    k=1;
    if row+k+1>size(edges_rep)
        break
    end
    while edges_rep(row+k,1)==edges_rep(row,1)
        k = k+1;
        if row+k>size(edges_rep)
            break
        end
    end
    [~,ord2]=sort(edges_rep(row:row+k-1,2),1);
    A = edges_rep(row:row+k-1,:);
    edges_rep(row:row+k-1,:) = A(ord2,:);
    row = row+k;
end

% construct the edges structure without repetitions
edges = [edges_rep(1,:) 0];
row = 1;
for i=2:size(edges_rep,1)
    if edges_rep(i,1:2)==edges(row,1:2)
        edges(row,4) = edges_rep(i,3);
        %edges(row,5) = edges(row,5)+edges_rep(i+1,5);
    else
        edges = [edges; [edges_rep(i,:) 0]];
        row = row+1;
    end
end

% compute edges midpoints
mid = [0.5*(nodes(edges(:,1),1)+nodes(edges(:,2),1)) 0.5*(nodes(edges(:,1),2)+nodes(edges(:,2),2))];


% add the edges' distances d_sig, length m_sig and the transmissivity d_sig/m_sig
% add also the distances d_{K,sig}, d_{L,sig}
edges = [edges zeros(size(edges,1),5)];
for e=1:size(edges,1)
    if edges(e,4)~=0
        d_sig = sqrt( (cc(edges(e,3),1)-cc(edges(e,4),1)).^2 + (cc(edges(e,3),2)-cc(edges(e,4),2)).^2 );
        dKsig = sqrt((mid(e,1)-cc(edges(e,3),1)).^2+...
            (mid(e,2)-cc(edges(e,3),2)).^2);
        dLsig = sqrt((mid(e,1)-cc(edges(e,4),1)).^2+...
            (mid(e,2)-cc(edges(e,4),2)).^2);
    else
        d_sig = sqrt( (cc(edges(e,3),1)-mid(e,1)).^2 + (cc(edges(e,3),2)-mid(e,2)).^2 );
        dKsig = sqrt((mid(e,1)-cc(edges(e,3),1)).^2+...
            (mid(e,2)-cc(edges(e,3),2)).^2);
        dLsig = 0;
    end
    m_sig = sqrt((nodes(edges(e,1),1)-nodes(edges(e,2),1)).^2+(nodes(edges(e,1),2)-nodes(edges(e,2),2)).^2);
    edges(e,5:9) = [d_sig m_sig m_sig/d_sig dKsig dLsig];
end



end