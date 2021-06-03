function [y] = mxv_jacobian(A,B1T,B2,C,ncellphi,x)
  n=size(B1T,1);
  m=size(B1T,2);

  y=zeros(n+m,1);

  v=B1T*x(n+1:n+m);
  
  y(1:n)     = A *x(1:n) + v;
  y(n+1:n+m) = B2*x(1:n) - C*x(n+1:n+m);

  imbalance_phi=0.0;
  imbalance_v=0.0;
  for i = 1:int64(n/ncellphi)
    imbalance_phi=max(imbalance_phi,sum(y((i-1)*ncellphi+1:i*ncellphi)));
    imbalance_v=max(imbalance_v,sum(v((i-1)*ncellphi+1:i*ncellphi)));
  end
  %fprintf('max imbalance A x* B1T*y = %1.4e, imb B1T*y = %1.4e \n',imbalance_phi,imbalance_v);
end
