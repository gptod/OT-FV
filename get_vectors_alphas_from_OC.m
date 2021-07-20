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

  if (~exist('approach','var') )
    approach=1;
  end
  
				% rhs
  f=-FOC.p;
  g=-FOC.r;
  h=-FOC.s;

  gtilde = g-JOC.rs*(JOC.ss\h);


  rho=-spdiags(JOC.ss,0);
  slack =-spdiags(JOC.sr,0);
  sizecell=-spdiags(JOC.rs(1:ncellrho,1:ncellrho),0);

  B2T=JOC.rp';
  B2=B2T';
  
  print=0;

  approach=controls.manipulation_approach;


  imbalances=zeros(Nt,1);
  betas=zeros(N,1);
  vi=zeros(Np,1);
  deltat=1/(N+1);
  for i=1:Nt
    imbalances(i)=sum(f(1+(i-1)*n:i*n));

    vi(:)=0;
    vi(1+(i-1)*n:i*n)=1.0;
    wi=B2*vi;
    if (0)
      fprintf('%d imb=%1.4e \n',i,imbalances(i))
      if (i<=N)
	fprintf('%d wi-sizecell/deltat=%1.4e \n',i,norm(wi(1+(i-1)*m:i*m)-sizecell/deltat))
      end
    end
    
  end
  betas(1)=imbalances(1);
  for i=2:N
    betas(i)= imbalances(i)+betas(i-1);
  end

  for i=1:N
    if (0)
      fprintf('%d beta=%1.4e\n',i,betas(i))
    end
  end

  

  if (approach==1)
    
				% allocations
    alphas=zeros(N);
    vectors=zeros(ncellphi,2*N);
    vectors_y=zeros(Nr,N);
    
    


    
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
				%sum(vectors(:,1+(i-1)*2))
				%sum(vectors(:,2+(i-1)*2))

				% vy=R*rho^i + M*slack^i
      vectors_y(:,i) = JOC.rr(:,1+(i-1)*m:i*m)*rho(1+(i-1)*m:i*m);
      vectors_y(1+(i-1)*m:i*m,i) = vectors_y(1+(i-1)*m:i*m,i)-(-sizecell.*slack(1+(i-1)*m:i*m));
      
      alphas(i)=(rho((1+(i-1)*m:i*m))'*g(1+(i-1)*m:i*m)) - sizecell'*h(1+(i-1)*m:i*m);
    end

  elseif(approach==2)
				% allocations
    alphas=zeros(N,1);
    vectors=zeros(ncellphi,2*N);
    imbalances=zeros(Nt,1);
    betas     =zeros(N,1);

    vectors_y=zeros(Nr,N);

    deltat=1/(N+1);

				% compute imbalances given by rhs
    for i=1:Nt
      imbalances(i)=sum(f(1+(i-1)*n:i*n));
      if (print)
	fprintf('%d imb=%1.4e \n',i,imbalances(i))
      end
    end
    betas(1)=imbalances(i)/deltat;
    for i=2:N
      betas(i)= imbalances(i)/deltat;
      betas(i) = betas(i)+betas(i-1);
      if (print)
	fprintf('%d beta=%1.4e\n',i,betas(i))
      end
    end


    if (1)
		  % M=diag(rho)*R+diag(s)*M
		  %rho=-spdiags(JOC.ss,0);
		  %slack =-spdiags(JOC.sr,0);
		  %sizecell=-spdiags(JOC.rs(1:ncellrho,1:ncellrho),0);
      M=JOC.ss*JOC.rr-(JOC.sr)*(JOC.rs);

      invM=sparse_inverse;

      ctrl_loc=ctrl_solver;
      index_agmg=1;
			      % passing index_agmg we will use agmg-i 
      ctrl_loc.init('direct',...
		    1e-12,...
		    100,...
		    0,...
		    0,...
		    [],...
		    index_agmg);
      invM.init(M,ctrl_loc);
    end

    vector_q=zeros(Nr,1);
				% g_i_tilda=diag(rho)*g - M * h
    g_tilde=rho.*g-JOC*h;

    inv_M_g_tilde=invM.apply(g_tilde);

    
    for i=1:N
			 % M(i,i)^{-1}=((R=0)*diag(rho_i)+m*diag(s_i))
%fprintf('(%1.2e <= S%d <= %1.2e) \n',min(slack(1+(i-1)*m:i*m)),i,max(slack(1+(i-1)*m:i*m)))
%invM=(1.0./(sizecell.*slack(1+(i-1)*m:i*m)));

				% support vector
      vector_q(:)=0;
      vector_q((i-1)*ncellrho+1:i*ncellrho)=1.0;

      wi=invM.apply(vector_qi);
      
				%
      vi=rho(1+(i-1)*ncellrho+1:i*ncellrho).*wi(1+(i-1)*ncellrho+1:i*ncellrho);
				% B2T(i,i)
      vectors(:,1+(i-1)*2)=B2T((i-1)*n+1:i*n,1+(i-1)*m:i*m)*vi;
				%sum(vectors(:,1+(i-1)*2))

				% B2T(i+1,i)
      vectors(:,2+(i-1)*2)=B2T(i*n+1:(i+1)*n,1+(i-1)*m:i*m)*vi;
				%sum(vectors(:,2+(i-1)*2))
      
			    % <p^i, M^{-1}(diag(rho) g - M h)>-beta(i)
      alphas(i)=vector_q'*inv_M_g_tilde-betas(i);
    end


    invM.kill();

  elseif(approach==3)
				% get B2 x -C y = ~g
    
    C = sparse(JOC.rr - JOC.rs*(JOC.ss\JOC.sr));
    %gtilde = g-JOC.rs*(JOC.ss\h);

	      % swap C sign for having standard saddle point notation 
    C = -C;

    CT=C';
    
				% allocations
    alphas=zeros(N,1);
    vectors=zeros(ncellphi,2*N);
    imbalances=zeros(Nt,1);
    betas     =zeros(N,1);

    vectors_y=zeros(Nr,N);

    deltat=1/(N+1);

				% compute imbalances given by rhs
   
    vi=ones(m,1);
    for i=1:N
				% B2T(i,i)*vi
      vectors(:,1+(i-1)*2)=B2T(1+(i-1)*n:i*n,1+(i-1)*m:i*m)*vi;
      
      
				% B2T(i+1,i)*vi
      vectors(:,2+(i-1)*2)=B2T(1+i*n:(i+1)*n,1+(i-1)*m:i*m)*vi;
      
				% q=-C(:,i)*wi
      vectors_y(:,i) = -CT(:,1+(i-1)*m:i*m)*vi;
      
      alphas(i)=vi'*(gtilde(1+(i-1)*m:i*m));
    end

  elseif(approach==4)
				% get B2 x -C y = ~g
    
    C = sparse(JOC.rr - JOC.rs*(JOC.ss\JOC.sr));

	      % swap C sign for having standard saddle point notation 
    C = -C;

    CT=C';
    
				% allocations
    alphas=zeros(N,1);
    vectors=zeros(ncellphi,2*N);

    vectors_y=zeros(Nr,N);

    deltat=1/(N+1);

    check_supports=0;
				% compute imbalances given by rhs

    kernels=zeros(Np,1);
    for i=1:Nt
      imbalances(i)=sum(f(1+(i-1)*n:i*n));
      if (check_supports)
	fprintf('%d imb=%1.4e \n',i,imbalances(i))
      end
      kernels(:)=0;
      kernels(1+(i-1)*n:i*n)=1.0;
      if (check_supports>1)
	%fprintf('%d norm(A*kernel)= %1.4e \n',i,norm(JOC.pp*kernels))
      
	support=blanks(N);
	B_times_kernel=B2*kernels;
	support(:)='_';
	for j=1:N
	  if (norm(B_times_kernel(1+(j-1)*m:j*m)-sizecell/deltat)<1e-13)
	    support(j)='+';
	  end
	  if (norm(B_times_kernel(1+(j-1)*m:j*m)+sizecell/deltat)<1e-13)
	    support(j)='-';
	  end
	end
	fprintf('B*kernels %03d %s \n',i,support);
      end
      
    end
    


    invCT=sparse_inverse;

    ctrl_loc=ctrl_solver;
    index_agmg=1;
			      % passing index_agmg we will use agmg-i 
    ctrl_loc.init('direct',...
		  1e-12,...
		  100,...
		  0,...
		  0,...
		  [],...
		  index_agmg);
    invCT.init(CT,ctrl_loc);
    
    
    vector=zeros(Nr,1);
    for i=1:N
      vector(:)=0;
      vector(1+(i-1)*m:i*m)=sizecell/deltat;
      wi=invCT.apply(vector);
      if (0)
	invCT.info_inverse.print();
      end
      B2T_times_wi=B2T*wi;

      if (check_supports)
	support=blanks(N+1);
	support(:)='_';
	for j=1:N+1
	  if (norm(B2T_times_wi(1+(j-1)*n:j*n))>1e-15)
	    if ( min(B2T_times_wi(1+(j-1)*n:j*n))>0.0)
	      support(j)='+';
	    elseif ( max(B2T_times_wi(1+(j-1)*n:j*n))<0.0)
	      support(j)='-';
	    else
	      support(j)='?';
	    end
	  end
	end
	fprintf('%03d %s \n',i,support);
      end
      
				% B2T(i,i)*vi
      vectors(:,1+(i-1)*2)=B2T_times_wi(1+(i-1)*n:i*n);
      %fprintf('%03d sum %1.2e norm %1.2e \n',i,sum(vectors(:,1+(i-1)*2)),norm(vectors(:,1+(i-1)*2)));
      
				% B2T(i+1,i)*vi
      vectors(:,2+(i-1)*2)=B2T_times_wi(1+i*n:(i+1)*n);
      %fprintf('%03d sum %1.2e norm %1.2e\n',i,sum(vectors(:,2+(i-1)*2)),norm(vectors(:,2+(i-1)*2)));
      
				% q=-C(:,i)*wi
      vectors_y(:,i) = -CT*wi;
      if (0)
	support=blanks(N);
	support(:)='_';
	for j=1:N
	  if (norm(vectors_y(1+(j-1)*m:j*m,i))>1e-15)
	    support(j)='o';
	  end
	end
	fprintf('%03d %s \n',i,support);
	fprintf('%03d norm(CT*wi+area)%1.2e \n',i,norm(vectors_y(:,i)+vector));
      end
      
      
      vectors_y(:,i)=0;

      alphas(i)=wi'*(gtilde)+betas(i);
    end

    invCT.kill();
    
  end


  if (controls.swap_sign)
    vectors=-vectors;
    vectors_y=-vectors_y;
    alphas=-alphas;
  end

  if (0)
    for i=1:N
      fprintf('(%1.4e <= p1 <= %1.4e) sum=%1.2e\n',min(vectors(:,1+(i-1)*2)),max(vectors(:,1+(i-1)*2)),sum(vectors(:,1+(i-1)*2)))
      fprintf('(%1.4e <= p2 <= %1.4e) sum=%1.2e \n',min(vectors(:,2+(i-1)*2)),max(vectors(:,2+(i-1)*2)),sum(vectors(:,2+(i-1)*2)))
    end

  end
end
