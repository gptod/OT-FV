function [nodes,cells] = get_mesh_2d(filename1,filename2,varargin)

% j1l and j1r number of junk left/right columns of the first file,
% j2l and j2r number of junk left/right columns of the second file,
% to be excluded from the mesh structure

if isempty(varargin)==0
    varargin = cell2mat(varargin);
    if length(varargin)==4
        j1l = varargin(1);
        j1r = varargin(2);
        j2l = varargin(3);
        j2r = varargin(4);
    else
        disp('error with the number of extra arguments in the function "get_mesh"')
    end
else
    j1l = 0;
    j1r = 0;
    j2l = 0;
    j2r = 0;
end
nodes = load(filename1);
nodes = nodes(:,j1l+1:end-j1r);

cells = read(filename2,j2l,j2r);
% add the number of vertices for each polygon
cells = [sum(cells~=0,2) cells];
% sort the vertices per each polygon
%for i=1:size(cells,1)
%    cells(i,2:1+cells(i,1)) = sort(cells(i,2:1+cells(i,1)));
%end
