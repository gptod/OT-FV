clear
close all

% construct the finer mesh from a coarse triangular mesh

name_f = 'subtri2_mesh1'; %name of the final structure

% load the coarse mesh
load('tri2_mesh1.mat');

if exist('cell_sig','var')==0
    % compute cells-edges structure
    maxn= max(cells(:,1));
    [cell_e,~,~] = str_cell_2d(maxn,ind,sigma,mid_edges,cc);
end

cells_f = zeros(3*ncell,5);
nodes_f = zeros(nnode+nsig+ncell,2);
nodes = [nodes zeros(nnode,1)];
mid_edges = [mid_edges zeros(nsig,1)];
cc = [cc zeros(ncell,1)];
int = 0;
I = sparse(3*ncell,ncell);

for k=1:ncell
    
    K = cells(k,2:4);
    E = cell_sig(k,2:4);
    
    for j=1:3
        
        
        ej = E(edges(E',1)==K(j) | edges(E',2)==K(j));
        
        if nodes(K(j),3)==0
            int = int+1;
            A = int;
            nodes(K(j),3) = int;
            nodes_f(int,:) = nodes(K(j),1:2);
        else
            A = nodes(K(j),3);
        end
        
        if mid_edges(ej(1),3)==0
            int = int+1;
            B = int;
            mid_edges(ej(1),3) = int;
            nodes_f(int,:) = mid_edges(ej(1),1:2);
        else
            B = mid_edges(ej(1),3);
        end
        
        if cc(k,3)==0
            int = int+1;
            C = int;
            cc(k,3) = int;
            nodes_f(int,:) = cc(k,1:2);
        else
            C = cc(k,3);
        end
        
        if mid_edges(ej(2),3)==0
            int = int+1;
            D = int;
            mid_edges(ej(2),3) = int;
            nodes_f(int,:) = mid_edges(ej(2),1:2);
        else
            D = mid_edges(ej(2),3);
        end
        
        cells_f(3*(k-1)+j,:) = [4 A B C D];
        
    end
    
    I(3*(k-1)+[1;2;3],k) = ones(3,1);
    
end

nodes(:,3) = [];
mid_edges(:,3) = [];
cc(:,3) = [];

nodes2h = nodes;
cells2h = cells;
edges2h = edges;
sigma2h = sigma;
ind2h = ind;
nodes = nodes_f;
cells = cells_f;
cc2h = cc;
mid_edges2h = mid_edges;
area2h = area;
ncell2h = ncell;
h2h = h;
clearvars -except nodes cells I name_f cc2h mid_edges2h area2h ncell2h nodes2h cells2h edges2h sigma2h ind2h h2h


nnode = size(nodes,1);
ncell = size(cells,1);
maxn= max(cells(:,1));


%% compute basic structures

% compute cells related properties
[cc,area,h] = mesh_2d(nodes,cells);

% compute edges structure
[sigma,edges,mid_edges] = str_sigma_2d(nodes,cells,cc);
nsig = size(sigma,1);

% compute indices
ind = indices(ncell,sigma);
nsig_in = length(ind.internal);


% %% compute final structures (optional)
% 
% % compute cells-edges structure
% [cell_sig,cell_sig_in,cell_dist] = str_cell_2d(maxn,ind,sigma,mid,cc);
% 
% % compute diamond cells structure
% [dnodes,dcells] = diamond_2d(nodes,cc,edges,sigma);

% print mesh
patch('Faces',edges(:,1:2),'Vertices',nodes)
axis square
axis off

%% saving

% save the mesh structure
save(name_f,'area','cc','cells','edges','sigma','h','ind','mid_edges','ncell','nsig','nsig_in','nnode','nodes','I',...
            'cc2h','mid_edges2h','area2h','ncell2h','nodes2h','cells2h','edges2h','sigma2h','ind2h','h2h')

