function [indeces] = set_grounding_node(A,ncellphi)
  Np=size(A,1);
  indeces=zeros(fix(Np/ncellphi),1);
  for i=1:fix(Np/ncellphi)
    [amax,imax]=max(abs(spdiags(A(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi),0)));
    indeces(i)=(i-1)*ncellphi+imax;
    %disp(indeces(i));
  end
end
