function [cc,area,h] = mesh(nodes,cells)

% INPUT:
% node: nodes' coordinate
% cell: topology
% OUTPUT:
% cc: circumcenters of the polygons
% area: area of the polygons
% h: meshsize of the mesh


% compute the circumcenters of the cells, namely the intersection of the
% segments orthogonal to the edges and passing through their midpoints
% take for this purpose two random edges...
vert = cells(:,2:end); % disregard the column specifying the number of vertices per cell
                      % attention: the vert mat contains possibly something
                      % different from vertices at the right end!!!
nvert = cells(:,1);

% compute the orthogonal vector to the first edge vertex_1-vertex_2
v = [nodes(vert(:,1),2)-nodes(vert(:,2),2) nodes(vert(:,2),1)-nodes(vert(:,1),1)];
% use the second edge vertex_2-vertex_3 if a node's x-ccordinate is equal to zero
corrv = find(v(:,1)==0);
v(corrv,:) = [nodes(vert(corrv,2),2)-nodes(vert(corrv,3),2) nodes(vert(corrv,3),1)-nodes(vert(corrv,2),1)];
if any(v(:,1)==0)
    error('there is a problem with the mesh')
end

% compute the orthogonal vector to the last edge vertex_1-vertex_nvert
w = zeros(size(vert,1),2);
for i=1:size(vert,1)
    w(i,:) = [nodes(vert(i,1),2)-nodes(vert(i,nvert(i)),2) nodes(vert(i,nvert(i)),1)-nodes(vert(i,1),1)];
end
% detect parallel edges (v//w) and change the edge w: use instead the edge
% vertex_(nvert-1)-vertex_nvert 
corrw = [];
for i=1:size(v)
    if w(i,:)-v(i,:)==[0 0]
        corrw = [corrw; i];
    end
end
for i=1:length(corrw)
    w(corrw(i),:) = [nodes(vert(corrw(i),nvert(corrw(i))-1),2)-nodes(vert(corrw(i),nvert(corrw(i))),2) nodes(vert(corrw(i),nvert(corrw(i))),1)-nodes(vert(corrw(i),nvert(corrw(i))-1),1)];
end

% compute the midpoints of the two edges considered
midv = [0.5*(nodes(vert(:,1),1)+nodes(vert(:,2),1)) 0.5*(nodes(vert(:,2),2)+nodes(vert(:,1),2))];
midv(corrv,:) = [0.5*(nodes(vert(corrv,2),1)+nodes(vert(corrv,3),1)) 0.5*(nodes(vert(corrv,3),2)+nodes(vert(corrv,2),2))];
midw = zeros(size(vert,1),2);
for i=1:size(vert,1)
    midw(i,:) = [0.5*(nodes(vert(i,1),1)+nodes(vert(i,nvert(i)),1)) 0.5*(nodes(vert(i,nvert(i)),2)+nodes(vert(i,1),2))];
end
for i=1:length(corrw)
    midw(corrw(i),:) = [0.5*(nodes(vert(corrw(i),nvert(corrw(i))-1),1)+nodes(vert(corrw(i),nvert(corrw(i))),1)) 0.5*(nodes(vert(corrw(i),nvert(corrw(i))),2)+nodes(vert(corrw(i),nvert(corrw(i))-1),2))];
end

% compute the point of intersection between the two vectors
beta = (midw(:,2)-midv(:,2)-v(:,2).*(midw(:,1)-midv(:,1))./v(:,1))./(v(:,2).*w(:,1)./v(:,1)-w(:,2));
cc = midw + beta.*w;


% a = [nodes(cells(:,1),1) nodes(cells(:,1),2) nodes(cells(:,1),3)];
% b = [nodes(cells(:,2),1) nodes(cells(:,2),2) nodes(cells(:,2),3)];
% c = [nodes(cells(:,3),1) nodes(cells(:,3),2) nodes(cells(:,3),3)];
% ac = c-a;
% ab = b-a;
% abxac = cross(ab,ac,2);
% abxacxab = cross(abxac,ab,2);
% acxabxac = cross(ac,abxac,2);
% cc = zeros(size(tri,1),3);
% for k=1:size(tri,1)
%     cc(k,:) = a(k,:) + (norm(ac(k,:),2)^2*abxacxab(k,:) + norm(ab(k,:),2)^2*acxabxac(k,:))/(2*norm(abxac(k,:),2)^2);
% end


% compute the area of the polygons
area = zeros(size(vert,1),1);
for i=1:size(vert,1)
    vert_i = vert(i,1:nvert(i));
    x_i = nodes(vert_i,:);
    area(i) = area_pol(nvert(i),x_i);
end

% for a manifold:
% for k=1:size(tri,1)
%     x1 = node(tri(k,1),1);
%     x2 = node(tri(k,2),1);
%     x3 = node(tri(k,3),1);
%     y1 = node(tri(k,1),2);
%     y2 = node(tri(k,2),2);
%     y3 = node(tri(k,3),2);
%     z1 = node(tri(k,1),3);
%     z2 = node(tri(k,2),3);
%     z3 = node(tri(k,3),3);
%     areak = 0.5*norm(cross([x2-x1; y2-y1; z2-z1],[x3-x1; y3-y1; z3-z1]),2);
%     area = [area; areak];
% end

% compute the meshsize: max(diam(polygon))
if any(nvert~=3)
    % to simplify the problem use the square root of the area instead of
    % the diameter 
    %h = max(sqrt(area));
    
    % for a convex polygon the diameter is equal to the maximum distance
    % between vertices
    h = 0;
    for i=1:size(vert,1)
        for j=1:nvert(i)
            for k=j+1:nvert(i)
                h = max(h, sqrt((nodes(vert(i,j),1)-nodes(vert(i,k),1)).^2+(nodes(vert(i,j),2)-nodes(vert(i,k),2)).^2) );
            end
        end
    end
else
    % if the polygons are triangle
    h1 = sqrt((nodes(vert(:,1),1)-nodes(vert(:,2),1)).^2+(nodes(vert(:,1),2)-nodes(vert(:,2),2)).^2);
    h2 = sqrt((nodes(vert(:,1),1)-nodes(vert(:,3),1)).^2+(nodes(vert(:,1),2)-nodes(vert(:,3),2)).^2);
    h3 = sqrt((nodes(vert(:,2),1)-nodes(vert(:,3),1)).^2+(nodes(vert(:,2),2)-nodes(vert(:,3),2)).^2);
    h = max(max([h1 h2 h3]));
end

end