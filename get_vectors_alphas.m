function [vectors,alphas] = get_vectors_alphas(sizecell,ncellphi,ncellrho,A,B1T,B2,C,invC,f,g)
  % find 2N vectors and N alphas  satisfing
  %
  % v1^t x_1 + v2^t x2 =a1
  % v3^t x_2 + v4^t x3 =a2
  % ....
  % v(2N-1)^t x_N + v(2N)^t x(N+1) =aN

  % from B2 x -C y = g we get 
  % y=C^{-1} (g- B2 x)
  %
  % We now use that and  (yi)' masses=0.0, thus
  % taking
  % vi=(0,..,masses, 0)
  %          ^
  %          |
  %          i-th block
  % we get 
  % <vi,y>=<vi,C^{-1} (g-B2x)>=0
  % alpha(i)=g^T vi
  % v1=B11 C^{-1} v1
  % v2=B12 C^{-1} v1
  % v3=B22 C^{-1} v2
  % v4=B23 C^{-1} v2

  % size
  N=fix(size(A,1)/ncellphi)-1;
  Nt=N+1;
  m=ncellrho;

  % allocations
  alphas=zeros(N);
  vectors(ncellphi,2*N);
  
  
  constant_masses=zerosas(g);
  for i=1:Nt
    constant_masses(1+(i-1)*m:i*m)=sizecell;
  end
  w=invC(constant_masses);

  for i=1:N
    vectors(:,1+(i-1)*2)=B2'(1+(i-1)*m:i*m,...
			     1+(i-1)*m:i*m)*...
			 w(1+(i-1)*m:i*m);
    vectors(:,2+(i-1)*2)=B2'(1+(i-1)*m:i*m,...
			     1+i*m:(i+1)*m)*...
			 w(1+(i-1)*m:i*m);
    alphas(i)=g(1+(i-1)*m:i*m)'*w(1+(i-1)*m:i*m);
  end

end
  
  
