function [ data ] = read_td_sequence( headbody,fid )
%UNTITnLED Summary of this function goes here
%   Detailed explanation goes here

if (headbody==1)
  tline = fgetl(fid);
  info=sscanf(tline,'%d %d  !*');
  dim=info(1)
  ndata=info(2)
  data=zeros(ndata,dim);
end

if (headbody==2)
  tline = fgetl(fid);
  nnz=int32(sscanf(tline,'%d!*'));
  
  for i=1:nnz;
    tline = fgetl(fid);
    info=sscanf(tline,'%d %f *');
    innz=int32(info(1));
    data(innz,1)=info(2);
  end
end

end


