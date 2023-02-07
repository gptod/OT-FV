function [output] = read(nomefile,varargin)

% jl and jr number of junk left/right columns of the file,
% to be excluded from the mesh structure

if isempty(varargin)==0
    varargin = cell2mat(varargin);
    if length(varargin)==2
        jl = varargin(1);
        jr = varargin(2);
    else
        disp('error with the number of extra arguments in the function "get_mesh"')
    end
else
    jl = 0;
    jr = 0;
end

fp=fopen(nomefile,'r');
output = [];
while ~feof(fp)
    str=fgetl(fp);
    data=sscanf(str,'%i'); data = data'; data = data(jl+1:end-jr);
    if isempty(data)==1
        % if a blank line is encountered end the process
        break
    elseif length(data)>size(output,2)
        % if a bigger polytope is detected, add blank columns to the
        % previous lines to fit data
        output = [output zeros(size(output,1),length(data)-size(output,2))];
    elseif length(data)<size(output,2)
        % if the number of actual vertices is less than tha maximum number
        % add zeros to the data to fit output
        data = [data zeros(1,size(output,2)-length(data))];
    end
    output = [output; data];
end
fclose(fp);


