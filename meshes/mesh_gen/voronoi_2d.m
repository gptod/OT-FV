function  [nodes_v,cells_v] = voronoi(nodes_t,edges_t,sigma_t,cc_t,mid,corners)

% Transform the Delaunay partitioning of a polygonal domain into a Voronoi tessellation
% The corners needs to be specified in order to be properly added

ncell_t = size(cc_t,1);
nnode_t = size(nodes_t,1);
nedge_t = size(edges_t,1);
nedgeb_t = size(sigma_t(sigma_t(:,2)==0),1);
nodes_v = [cc_t;zeros(nedgeb_t,2);nodes_t(corners,:)];

Nmax = 4;
ncells_v = nnode_t;
cells_v = zeros(ncells_v,Nmax+1);
count = 0;

for e=1:nedge_t
    
    A = edges_t(e,1);
    B = edges_t(e,2);
    K = sigma_t(e,1);
    L = sigma_t(e,2);
    
    if L~=0
        
        cells_v(A,1) = cells_v(A,1)+2;
        if cells_v(A,1) > Nmax
            cells_v = [cells_v zeros(nnode_t,cells_v(A,1)-Nmax)];
            Nmax = cells_v(A,1);
        end
        cells_v(A,cells_v(A,1)) = K;
        cells_v(A,1+cells_v(A,1)) = L;
        
        cells_v(B,1) = cells_v(B,1)+2;
        if cells_v(B,1) > Nmax
            cells_v = [cells_v zeros(nnode_t,cells_v(B,1)-Nmax)];
            Nmax = cells_v(B,1);
        end
        cells_v(B,cells_v(B,1)) = K;
        cells_v(B,1+cells_v(B,1)) = L;
        
    else
        
        cells_v(A,1) = cells_v(A,1)+2;
        if cells_v(A,1) > Nmax
            cells_v = [cells_v zeros(nnode_t,cells_v(A,1)-Nmax)];
            Nmax = cells_v(A,1);
        end
        cells_v(A,cells_v(A,1)) = K;
        count = count+1;
        nodes_v(ncell_t+count,:) = mid(e,:);
        cells_v(A,1+cells_v(A,1)) = ncell_t+count;
        
        cells_v(B,1) = cells_v(B,1)+2;
        if cells_v(B,1) > Nmax
            cells_v = [cells_v zeros(nnode_t,cells_v(B,1)-Nmax)];
            Nmax = cells_v(B,1);
        end
        cells_v(B,cells_v(B,1)) = K;
        cells_v(B,1+cells_v(B,1)) = ncell_t+count;
    
    end
    
    % add the corners
    if any(A==corners)
        cells_v(A,1) = cells_v(A,1)+1;
        if cells_v(A,1) > Nmax
            cells_v = [cells_v zeros(nnode_t,cells_v(A,1)-Nmax)];
            Nmax = cells_v(A,1);
        end
        cells_v(A,1+cells_v(A,1)) = ncell_t+nedgeb_t+find(A==corners);
    elseif any(B==corners)
        cells_v(B,1) = cells_v(B,1)+1;
        if cells_v(B,1) > Nmax
            cells_v = [cells_v zeros(nnode_t,cells_v(B,1)-Nmax)];
            Nmax = cells_v(B,1);
        end
        cells_v(B,1+cells_v(B,1)) = ncell_t+nedgeb_t+find(B==corners);
    end
        
        
end

% eliminate doubles
cells_v2 = zeros(size(cells_v));
for k=1:nnode_t
    [~,I,~] = unique(cells_v(k,2:1+cells_v(k,1)),'stable');
    cells_v2(k,1:length(I)+1) = [length(I) cells_v(k,1+I)];
end
cells_v = cells_v2;
maxI = max(cells_v(:,1));
cells_v = cells_v(:,1:maxI+1);


% sort the nodes of the cells
for k=1:ncells_v
    
    % The voronoi tessellation is made of convex cells
    % Compute the convex hull of the nodes of the cell thanks to matlab
    % function convhull
    % The convex hull K is expressed in terms of a vector of point indices 
    % arranged in a counter-clockwise cycle around the hull.
    N = cells_v(k,2:cells_v(k,1)+1);
    K = convhull(nodes_v(N,:));
    cells_v(k,2:cells_v(k,1)+1) = N(K(1:end-1));
    
end

% replace zeros with the first element
for k=1:nnode_t
    n = cells_v(k,1);
    cells_v(k,n+2:end) = cells_v(k,2);
end






