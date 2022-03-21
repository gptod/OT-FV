function x_new = shift_solution(x,vectors,alphas)
  % we shift (x2,...,xN) but not x1
  % x=(x1; x2; ...x(N+1)) so that
	% v1^t x^1 = 0.0
  % v1^t x_1 + v2^t x2 =a1
  % v3^t x_2 + v4^t x3 =a2
  % ....
  % v(2N-1)^t x_N + v(2N)^t x(N+1) =aN
  %
  % It is a sort of forward propagation of the correction

  % arrays size
  ncellphi=size(vectors,1);
  N=size(alphas,1);
  %N=size(vectors,2);

  
  for i=1:N
    inext=i+1
    % for the first 
    % alpha=a1-v1^t x_1
    alpha=alphas(i)-vectors(1+(i-1)*2)'*x(1+(i-1)*ncellphi:i*ncellphi);
    % v2^x2+shift sum(v2)=alpha
    v=vectors(:,2+(i-1)*2);
    sumv=sum(v);
    shift=(alpha-v'*x(1+i*ncellphi:(i+1)*ncellphi))/sumv
    % shift by constant
    x_new(1+i*ncellphi:(i+1)*ncellphi)=x(1+i*ncellphi:(i+1)*ncellphi)+shift;
  end
end
  
