function [vectors,alphas] = get_vectors_alphas(sizecell,ncellphi,ncellrho,A,B1T,B2,C,invC,f,g)
  % find 2N vectors and N alphas  satisfing
  %
  % v1^t x_1 + v2^t x2 =a1
  % v3^t x_2 + v4^t x3 =a2
  % ....
  % v(2N-1)^t x_N + v(2N)^t x(N+1) =aN

  % from B2 x -C y = g we get 
  % y=C^{-1} (B2 x-g)
  %
  % We now use that and  (yi)' masses=b_i=sum(f^i)*deltat, thus
  % taking
  % vi=(0,..,masses, 0)
  %          ^
  %          |
  %          i-th block
  % we get 
  % b_i=<vi,y>=<vi,C^{-1} (g-B2x)>=0
  % alpha(i)=<C^{-T} vi,g>-bi
  % v1=B11 C^{-1} v1
  % v2=B12 C^{-1} v1
  % v3=B22 C^{-1} v2
  % v4=B23 C^{-1} v2

  % size
  N=size(A,1)/ncellphi-1
  Nt=N+1;
  ncellrho;
  n=ncellphi;
  m=ncellrho;
  deltat=1.0/Nt;

  % allocations
  alphas=zeros(N,1);
  vectors=zeros(ncellphi,2*N);
  imbalances=zeros(Nt,1);
  betas     =zeros(N,1);



  % compute imbalances given by rhs
  for i=1:Nt
    imbalances(i)=sum(f(1+(i-1)*n:i*n));
    fprintf('%d imb=%1.4e \n',i,imbalances(i))
  end
  betas(1)=imbalances(i)/deltat;
  for i=2:N
    betas(i)= imbalances(i)/deltat;
    betas(i) = betas(i)+betas(i-1);
    fprintf('%d beta=%1.4e\n',i,betas(i))
  end

  
  
  % to optimize
  disp('B2T')
  B2T=B2';



  pvect=zeros(size(f,1),1);
  constant_masses=zeros(size(g,1),1);
  for i=1:N
    %fprintf('w%d \n',i)
    % da ottimizzare
    constant_masses(:)=0.0;
    constant_masses(1+(i-1)*m:i*m)=sizecell;
    w=invC(constant_masses);
    % check support
    %fprintf('C^{-1}w%d \n',i)
    %for j=1:N
    %  fprintf('%d %f\n',j,norm(w(1+(j-1)*m:j*m)))
    %end
    pvect=B2T*w;

    % check support
    %fprintf('B2 TC^{-1}w%d \n',i)
    %for j=1:N+1
    %  fprintf('%d %f\n',j,norm(pvect(1+(j-1)*n:j*n)));
    %end
    vectors(:,1+(i-1)*2)=pvect(1+(i-1)*n:i    *n);
    vectors(:,2+(i-1)*2)=pvect(1+i    *n:(i+1)*n);
    alphas(i)=g'*w-betas(i);
  end

  
  
end
  
  
