function [y] = mxv_jacobian(x,A,B1T,B2,C,ncellphi,kernel,rhs)
  n=size(B1T,1);
  m=size(B1T,2);

  y=zeros(n+m,1);

  v=B1T*x(n+1:n+m);


  
  verbose=0;
  if (verbose)
    print_imbalance(x(1:n),ncellphi)
  end
  
  %y(1:n)     = A *projector(x(1:n),kernel) + v;
  y(1:n)     = A *x(1:n) + v;
  y(n+1:n+m) = B2*x(1:n) - C*x(n+1:n+m);

  imbalance_phi=0.0;
  imbalance_v=0.0;

  if (verbose)
    print_imbalance(y(1:n),ncellphi)
  end

  if (verbose)
    if (exist('rhs','var') )
      fprintf(' res=%1.4e \n', norm(y-rhs)/norm(rhs))
      print_imbalance(rhs(1:n),ncellphi)
    end
  end
end
