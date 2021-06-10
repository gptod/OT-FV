function [ data ] = read_td( file )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fid = fopen(file);

tline = fgetl(fid);
info=sscanf(tline,'%d %d  !*');
dim=info(1);
ndata=info(2);
data=zeros(ndata,dim);


tline = fgetl(fid);
tline = fgetl(fid);
nnz=int32(sscanf(tline,'%d!*'));

for i=1:nnz;
    tline = fgetl(fid);
    info=sscanf(tline,'%d %f *');
    innz=int32(info(1));
    data(innz)=info(2);
end

fclose(fid);

end


