function [dense_S]=form_minus_SchurCA(inv_A,B1T,B2,C)
  Nr=size(C,1);
  dense_S=zeros(Nr,Nr);
  fprintf('building ~S=C+B1T invA B2')
  fprintf('%4.2f',0)
  for i=1:Nr
    if (mod(i,100)==0)
      fprintf('\b\b\b\b\b\b%4.2f ',i/Nr*100)
    end
    dense_S(:,i)=B2*(inv_A(full(B1T(:,i))));
  end
  
  dense_S=dense_S+C;
end
