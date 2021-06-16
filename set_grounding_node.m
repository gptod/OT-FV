function [indeces_global,indeces_local] = set_grounding_node(A,ncellphi)
  Np=size(A,1);
  indeces_global=zeros(fix(Np/ncellphi),1);
  indeces_global=zeros(fix(Np/ncellphi),1);
  for i=1:fix(Np/ncellphi)
    [amax,imax]=max(abs(spdiags(A(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi),0)));
    indeces_local(i)=imax;
    indeces_global(i)=(i-1)*ncellphi+imax;
    
  end
end
