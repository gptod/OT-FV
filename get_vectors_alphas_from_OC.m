function [vectors,vectors_y,alphas] = get_vectors_alphas_from_OC(JOC,FOC,controls)
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

  ncellphi=JOC.ncellphi;
  ncellrho=JOC.ncellrho;
  
  N=JOC.ntimestep;
  Nt=N+1;
  n=ncellphi;
  m=ncellrho;
  Np=Nt*ncellphi;
  Nr=N*ncellrho;

  
  
  % rhs
  f=-FOC.p;
  g=-FOC.r;
  h=-FOC.s;

  
 
  % allocations
  alphas=zeros(N);
  vectors=zeros(ncellphi,2*N);
  vectors_y=zeros(Nr,N);
  
  

  rho=-spdiags(JOC.ss,0);
  slack =-spdiags(JOC.sr,0);
  sizecell=-spdiags(JOC.rs(1:ncellrho,1:ncellrho),0);

  B2T=JOC.rp';
  
  for i=1:N
    % c1=B11^T rho^1
    % c2=B12^T rho^2
    % B2T*diag(rho^i)*1
    %fprintf('%d nrow %d,%d<=%d ncol %d,%d<=%d\n', i,1+(i-1)*n,i*n,Np,1+(i-1)*m,i*m,Nr)
    vectors(:,1+(i-1)*2)=B2T(1+(i-1)*n:i*n,...
			     1+(i-1)*m:i*m)*...
			 rho(1+(i-1)*m:i*m);
    %fprintf('%d nrow %d,%d<=%d ncol %d,%d<=%d\n', i,1+i*n,(i+1)*n,Np,1+(i-1)*m,i*m,Nr)
    vectors(:,2+(i-1)*2)=B2T(1+i*n:(i+1)*n,...
			     1+(i-1)*m:i*m)*...
			 rho(1+(i-1)*m:i*m);
    % D(rho) * (B2 x + Ry - M z) =D(rho) g
    % D(m)(-D(s) y-D(rho)z)=D(m)*h

    % vy=R*rho^i + M*slack^i
    vectors_y(:,i) = JOC.rr(:,1+(i-1)*m:i*m)*rho(1+(i-1)*m:i*m);
    vectors_y(1+(i-1)*m:i*m,i) = vectors_y(1+(i-1)*m:i*m,i)-(-sizecell.*slack(1+(i-1)*m:i*m));
    
    alphas(i)=(rho((1+(i-1)*m:i*m))'*g(1+(i-1)*m:i*m)) - sizecell'*h(1+(i-1)*m:i*m);
  end

end
  
  
