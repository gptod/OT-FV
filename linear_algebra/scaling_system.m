function [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs)

  Np=size(A,1);
  Nr=size(C,1);
  
  diagA_scaling = sparse(1:Np,1:Np,(1.0./sqrt(spdiags(A,0)))',Np,Np);
  diagC_scaling=sparse(1:Nr,1:Nr,(1.0./sqrt(spdiags(C,0)))',Nr,Nr);

  if (0)
    fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)));
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  end
  
  A=  diagA_scaling*A*  diagA_scaling;
  B1T=diagA_scaling*B1T*diagC_scaling;
  B2= diagC_scaling*B2* diagA_scaling;
  C=  diagC_scaling*C*  diagC_scaling;

%fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
%fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)))

  rhs(1:Np)      =diagA_scaling*rhs(1:Np);
  rhs(1+Np:Np+Nr)=diagC_scaling*rhs(1+Np:Np+Nr);

end
