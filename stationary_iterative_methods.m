function [x,k,res,resv] = stationary_iterative_methods(operator,rhs,x0,eps,kmax,prec)
%matrix version
  n=size(rhs,1);
k=0;
x=x0;
nb=norm(rhs);
residuum=-(operator(x)-rhs);
res=norm(residuum)/nb
resv=res;
while res>eps && k<kmax
  k=k+1;
  p=prec(residuum);
  x=x+p;
  residuum=-(operator(x)-rhs);
  res=norm(residuum)/nb;
  fprintf('%d %9.2e %9.2e \n',k,norm(p),res)
  resv = [resv;res];
end
