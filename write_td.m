function [] = write_td( filename ,data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ndata=size(data,1);
fileID = fopen(filename,'w');
fprintf(fileID,'1 %d \n',ndata);
fprintf(fileID,'time 0.0\n');
fprintf(fileID,'%d \n',ndata);
for i=1:ndata
    fprintf(fileID,'%d %16.8e \n',i, data(i));
end
fprintf(fileID,'time 1e30 \n');
fclose(fileID);
end


