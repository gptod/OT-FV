function [A,B1T,B2,C,f1,f2,W_mat]=scale_saddle_point(A,B1T,B2,C,f1,f2,W_mat,diagA_scaling,diagC_scaling)

  Np=size(A,1);
  Nr=size(C,1);
    
  A=  diagA_scaling*A*  diagA_scaling;
  B1T=diagA_scaling*B1T*diagC_scaling;
  B2= diagC_scaling*B2* diagA_scaling;
  C=  diagC_scaling*C*  diagC_scaling;

%fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
%fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)))

  f1 = diagA_scaling * f1;
  f2 = diagC_scaling * f2;

	W_mat = W_mat * diagC_scaling;

end
