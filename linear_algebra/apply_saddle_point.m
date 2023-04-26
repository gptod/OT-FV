function [y] = apply_saddle_point(x,A,B1T,B2,C,n,m,N)

  y=zeros(n+m,1);
	y(1:n)     = A(x(1:n))  + B1T(x(n+1:n+m));
	%print_imbalance(y(1:n),n/(N+1))
  y(n+1:n+m) = B2(x(1:n)) - C(x(n+1:n+m));

end
