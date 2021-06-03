function [out] = projector(v,kernel)

  dim=size(kernel,1);

  
  out=v;
				%maxs=0.0;
  for i=1:fix(size(v,1)/dim)
    scalin=out(1+(i-1)*dim:i*dim)'*kernel;
    %fprintf('%d scalin = 1.4e \n',i,scalin);
    %maxs=max(maxs,scalin);
    out(1+(i-1)*dim:i*dim)=...
    out(1+(i-1)*dim:i*dim) - scalin/(kernel'*kernel) * kernel;
    scalout=out(1+(i-1)*dim:i*dim)'*kernel;
  end
end
