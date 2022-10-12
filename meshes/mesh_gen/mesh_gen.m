clear
close all

% generate 2D mesh structures starting from coordinates and topology
% INPUT: 2D coordinates and topology
%        node structure: [x-coordinate, y-coordinate]
%        topology structure: [vertex 1, vertex 2, ...]
% OUTPU: mesh structure

fn1 = 'tri2/mesh1/coord.txt'; %coordinates
fn2 = 'tri2/mesh1/topol.txt'; %topology
name = 'tri2_mesh1'; %name of the final structure

%% load coordinates and topology
% use the last four input variables to erase useless columns from the two
% files...
%[nodes,cells] = get_mesh(fn1,fn2,1,0,1,1);
[nodes,cells] = get_mesh(fn1,fn2,0,0,1,0);
%[nodes,cells] = get_mesh(fn1,fn2,0,0,0,0);
nnode = size(nodes,1);
ncell = size(cells,1);
maxn= max(cells(:,1));

% shift coordinates (if needed)
% nodes(:,1) = (nodes(:,1)+1)*0.5;
% nodes(:,2) = (nodes(:,2)+1)*0.5;
%nodes(:,1) = (nodes(:,1))*1.5;
%nodes(:,2) = (nodes(:,2))*1.5;
%nodes(:,1) = nodes(:,1)*2;
%nodes(:,2) = nodes(:,2)*2;


%% compute basic structures

% compute cells related properties
[cc,area,h] = mesh(nodes,cells);

% compute edges structure
[edges,mid] = str_edge(nodes,cells,cc);
nedge = size(edges,1);

% compute indices
ind = indices(ncell,edges);
nei = length(ind.internal);


% %% reorder the elements in order to reduce the bandwidth (optional)
% 
% % compute the stiffness matrix
% A = FVstiff(ind,edges,ncell);
% % compute the permutation that reduces the bandwidth
% p = symrcm(A); p =p';
% % reorder the cells 
% cells = cells(p,:);
% %clearvars -except nodes cells nnode ncell maxn nedge nei name
% % recompute the structures
% cc = cc(p,:); area = area(p);
% [edges,mid] = str_edge(nodes,cells,cc);
% ind = indices(ncell,edges);


%% compute final structures (optional)

% compute cells-edges structure
%[cell_e,~,~] = str_cell(maxn,ind,edges,mid,cc);

% % compute diamond cells structure
% [dnodes,dcells] = diamond(nnode,nedge,nodes,cc,edges);


%% saving

% save the mesh structure
save(name,'area','cc','cells','edges','h','ind','mid','ncell','nedge','nei','nnode','nodes')


% %% check validity of the discretization for triangulations
% 
% %check for obtuse angles
% obtuse=[];
% bar = 1/3*(nodes(cells(:,2),:) + nodes(cells(:,3),:) + nodes(cells(:,4),:));
% 
% for k=1:ncell
%     
%     xe1 = mid(cell_e(k,2),:);
%     xe2 = mid(cell_e(k,3),:);
%     xe3 = mid(cell_e(k,4),:);
%     
%     xk = cc(k,:);
%     xk_grav=bar(k,:);
%     
%     scal=zeros(1,3);
%     
%     scal(1) = (xk-xe1)*(xk_grav-xe1)';
%     scal(2) = (xk-xe2)*(xk_grav-xe2)';
%     scal(3) = (xk-xe3)*(xk_grav-xe3)';
%     
%     if min(scal,[],2)<0
%         obtuse = [obtuse,k];
%     end
%     
% end
% 
% if isempty(obtuse)==0
%     disp('obtuse angles!!!')
% end

% print mesh
patch('Faces',edges(:,1:2),'Vertices',nodes)
axis square
axis off

