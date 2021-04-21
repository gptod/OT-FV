function [d, sol_stat, normd] = solvesys(JF,F,indc,sol,ctrl_inner,ctrl_outer,compute_eigen,logID)

% solve the linear system d = -JF\F, reducing the system by eliminating the
% additional variable s and imposing the condition d(indc) = 0 to fix the 
% constant of the potential dp


Np = length(F.p);
Nr = length(F.r);
Ns = length(F.s);
N = JF.ntimestep;
Nt=N+1;
ncellphi=JF.ncellphi;
ncellrho=JF.ncellrho;
inner_nequ=0;

flag=0;
iter=0;
res=1e30;


% rhs
f=-F.p;
g=-F.r;
h=-F.s;

% Solve with the Schur complement A|C, using Matlab's backslash (or agmg)
swap_sign=1;

				% reduced system to phi rho variables
A = sparse(JF.pp); B1T = sparse(JF.pr); B2 = sparse(JF.rp);
C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
f1 = f;
f2 = g-JF.rs*(JF.ss\h);

% swap C sign for having standard saddle point notation 
C = -C;
if (swap_sign)
  A=-A;
  B1T=-B1T;
  B2=-B2;
  C=-C;
  f1=-f1;
  f2=-f2;
end

% assembly full system
jacobian = [A B1T; B2 -C];

% assembly rhs
rhs=[f1;f2];


if sol==1
    
    % Solve using Matlab's backslash
    
    %eps = 1e-8;
    %eId = eps*speye(size(JF.pp));
    Sys = [JF.pp JF.pr; JF.rp JF.rr-JF.rs*(JF.ss\JF.sr)];
    Sys(indc,:) = sparse(1,Np+Nr); Sys(:,indc) = sparse(Np+Nr,1); Sys(indc,indc) = 1;
    f1 = F.p; f1(indc) = 0;
    f2 = F.r-JF.rs*(JF.ss\F.s);
    d = Sys\[-f1;-f2];
    dp = d(1:Np); dr = d(Np+1:end);
    ds = JF.ss\(-F.s-JF.sr*dr);
    d = [dp; dr; ds];
    
    normd = norm(d);
    
    sol_stat = struct('flag',0,'relres',[],'iter',[]);
    
elseif sol==2
    
    % Solve with the Schur complement A|C, using Matlab's backslash (or agmg)
  bs=tic;
    A = JF.pp; B1 = JF.pr; B2 = JF.rp; invC = -((JF.sr*JF.rs)\JF.ss);
    f1 = F.p;
    f2 = F.r-JF.rs*(JF.ss\F.s);
    S = -A+B1*(invC*B2); S = S(2:end,2:end);
    fp = f1-B1*(invC*f2); fp = fp(2:end);
    dp = S\fp;
    outer_cpu=toc(bs);
    outer_iter=0;
    sol_stat = struct('flag',0,'relres',[],'iter',[]);
    
    %[dp,flag,res,iter] = agmg(S,fp,0,1e-6,100,0,[],0);
    %sol_stat = struct('flag',flag,'relres',res,'iter',iter);

    inner_iter=0
    inner_cpu=0
    
    dp = [0; dp];
    dr = invC*(-f2-B2*dp);
    ds = JF.ss\(-F.s-JF.sr*dr);
    d = [dp; dr; ds];
    
    normd = norm(d);
    
elseif sol==3
    
    eta = 1e-5;
    
    % Solve the system with the flexible GMRES and the left (right) block
    % preconditioner P given by P=[S 0; B Cp] (P=[S B'; 0 Cp]),
    % with S=A-B^T*Cp^{-1}*B and Cp=diag(C)
    
    C = JF.rr-JF.rs*(JF.ss\JF.sr);
    Sys = [JF.pp JF.pr; JF.rp C];
    Sys(indc,:) = sparse(1,Np+Nr); Sys(:,indc) = sparse(Np+Nr,1); Sys(indc,indc) = 1;
    b = [-F.p; -F.r+JF.rs*(JF.ss\F.s)]; b(indc) = 0;
    
    % % scaling
    % diagSys = abs(diag(Sys)).^(-0.5); diagSys(diagSys<1e-6) = 1;
    % diagSys = spdiags(diagSys,0,Np+Nr,Np+Nr);
    % Sys = diagSys*Sys*diagSys; b = diagSys*b;
    
    Ap = Sys(1:Np,1:Np);
    B1p = Sys(1:Np,Np+1:Np+Nr);
    B2p = Sys(Np+1:Np+Nr,1:Np); 
    Cp = spdiags(diag(C),0,Nr,Nr);
    
    % Solve the system
    [d,iter] = fgmres(Sys,b,eta,'max_iters',1,'restart',ceil(sqrt(size(Sys,1))),...
        'verb',0,'P',@(x,tol) blockP(x,tol,Ap,B1p,B2p,Cp));
        
    %d = diagSys*d;
    
    dp = d(1:Np); dr = d(Np+1:end);
    ds = JF.ss\(-F.s-JF.sr*dr);
    d = [dp; dr; ds];
    
    normd = norm(d);
    
    sol_stat = struct('flag',0,'relres',[],'iter',iter);
    
elseif sol==4
      
    % Solve the augmented system of equations with the flexible GMRES and the
    % left (right) block preconditioner P given by P=[S 0; B2 CS]
    %(P=[S B1; 0 Cp]), with S=A-B1*Cp^{-1}*B2 and Cp=diag(C)
    
    C = JF.rr-JF.rs*(JF.ss\JF.sr);
    gamma = 1;
    A_aug = JF.pp+gamma*JF.pr*JF.rp; B1_aug = JF.pr-gamma*JF.pr*C;
    Sys = [A_aug B1_aug; JF.rp C];
    Sys(indc,:) = sparse(1,Np+Nr); Sys(:,indc) = sparse(Np+Nr,1); Sys(indc,indc) = 1;
    f1 = -F.p;
    f2 = -F.r+JF.rs*(JF.ss\F.s);
    b = [f1+gamma*JF.pr*f2; f2]; b(indc) = 0;
    
    Ap = Sys(1:Np,1:Np);
    B1p = Sys(1:Np,Np+1:Np+Nr);
    B2p = Sys(Np+1:Np+Nr,1:Np); 
    Cp = spdiags(diag(C),0,Nr,Nr);
    
    [d,iter] = fgmres(Sys,b,1e-3,'max_iters',size(Sys,1),...
                                  'restart',size(Sys,1),'verb',0,'P',@(x,tol) blockP(x,tol,Ap,B1p,B2p,Cp));
    flag = 0;
    relres = 0;
    
    dp = d(1:Np); dr = d(Np+1:end);
    ds = JF.ss\(-F.s-JF.sr*dr);
    d = [dp; dr; ds];
    
    normd = norm(d);
    
    sol_stat = struct('flag',flag,'relres',relres,'iter',iter);
    
    
elseif sol==5
    
    % Solve the system with the GMRES and the sparse ilu decomposition as
    % preconditioner
    
    Sys = [JF.pp JF.pr; JF.rp JF.rr-JF.rs*(JF.ss\JF.sr)];
    Sys(indc,:) = sparse(1,Np+Nr); Sys(:,indc) = sparse(Np+Nr,1); Sys(indc,indc) = 1;
    b = [-F.p; -F.r+JF.rs*(JF.ss\F.s)]; b(indc) = 0;
    
    [L,U] = ilu(Sys);
    
    [d,flag,relres,iter] = gmres(Sys,b,[],1e-3,size(Sys,1),L,U);
    
    dp = d(1:Np); dr = d(Np+1:end);
    ds = JF.ss\(-F.s-JF.sr*dr);
    d = [dp; dr; ds];
    
    normd = norm(d);
    
    sol_stat = struct('flag',flag,'relres',relres,'iter',iter);

elseif sol==6
  % rhs
  f=-F.p;
  g=-F.r;
  h=-F.s;
  
  % Solve with the Schur complement A|C, using Matlab's backslash (or agmg)
  swap_sign=1;

  % reduced system to phi rho variables
  A = JF.pp; B1T = JF.pr; B2 = JF.rp;
  C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
  f1 = f;
  f2 = g-JF.rs*(JF.ss\h);
  
  % swap C sign for having standard saddle point notation 
  C = -C;
  if (swap_sign)
    A=-A;
    B1T=-B1T;
    B2=-B2;
    C=-C;
    f1=-f1;
    f2=-f2;
  end

  %J=[A, B1T;  B2, -C];
  %norm(J*d(1:Np+Nr)-[f1;f2])

  % start solution
  % invert C
  fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  invC = sparse(1:Nr,1:Nr,(1.0./spdiags(C))',Nr,Nr);
 
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invC*B2);
  fp = f1+B1T*(invC*f2);

  %norm(S*d(1:Np)-fp)
  
  % ground the solution
  S(indc,:) = sparse(1,Np);
  S(:,indc) = sparse(Np,1);
  S(indc,indc) = 1;
  fp(indc)     = 0;
  
    
  % solve with multigrid
  %[dp,flag,res,iter] = agmg(S,fp,0,ctrl.tolerance,ctrl.itermax,0,[],0);
  invS=sparse_inverse;
  invS.init(S,ctrl_inner);
  dp=invS.apply(fp);
  %invS.info.print();
  flag=invS.info.flag;
  relres=invS.info.res;
  iter=invS.info.iter;
  invS.kill();

  %res;norm(S*dp-fp)

  % get solution of system JF d  = F  
  dr = invC*(B2*dp-f2);

  %J=[A, B1T;  B2, -C];
  %norm(J*[dp;dr]-[f1;f2])

  ds = JF.ss\(h-JF.sr*dr);
  d = [dp; dr; ds];
    
  normd = norm(d);

elseif sol==7
  % rhs
  f=-F.p;
  g=-F.r;
  h=-F.s;
  
  % Solve with the Schur complement A|C, using Matlab's backslash (or agmg)
  swap_sign=1;

  % reduced system to phi rho variables
  A = JF.pp; B1T = JF.pr; B2 = JF.rp;
  C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
  f1 = f;
  f2 = g-JF.rs*(JF.ss\h);
  
  % swap C sign for having standard saddle point notation 
  C = -C;
  if (swap_sign)
    A=-A;
    B1T=-B1T;
    B2=-B2;
    C=-C;
    f1=-f1;
    f2=-f2;
  end

  %S=compute_SchurCA(A,B1T,B2,C,Nt,ncellrho,ncellphi)
  %S=A+B1T*(C\B2);
  
  %J=[A, B1T;  B2, -C];
  %norm(J*d(1:Np+Nr)-[f1;f2])

  % start solution
  % invert C
  diagC=spdiags(C);
  %fprintf('(%9.4e <= C <= %18.12e) \n',min(diagC),max(diagC))
  invdiagC = sparse(1:Nr,1:Nr,(1.0./diagC)',Nr,Nr);
 
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invdiagC*B2);
  fp = f1+B1T*(invdiagC*f2);

 

  %norm(S*d(1:Np)-fp)
  
  % ground the solution
  S(indc,:) = sparse(1,Np);
  S(:,indc) = sparse(Np,1);
  S(indc,indc) = 1;
  fp(indc)     = 0;

  if ( compute_eigen )
    disp("START EIG")
    eigenvalues=eig(full(S));
    disp("END EIG")
    for i=1:4
      fprintf('(%9.4e <= eigvalues <= %18.12e) \n',eigenvalues(i),eigenvalues(Np-i+1))
    end
  end
    
  % solve with multigrid
  %ctrl.info();
  invS=sparse_inverse;
  invS.init(S,ctrl_inner);
  dp=invS.apply(fp);
  invS.info.print();
  flag=invS.info.flag;
  relres=invS.info.res;
  iter=invS.info.iter;
  invS.kill();
  
  %res;norm(S*dp-fp)

  % get solution of system JF d  = F  
  dr = invdiagC*(B2*dp-f2);

  %J=[A, B1T;  B2, -C];
  %norm(J*[dp;dr]-[f1;f2])

  ds = JF.ss\(h-JF.sr*dr);
  d = [dp; dr; ds];
    
  normd = norm(d);

elseif sol==8
  % Invert with respect to C block and try to "invert" S=A+BT*C^{-1} B
 
  % start solution
  % invert C
  fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  invC = sparse(1:Nr,1:Nr,(1.0./spdiags(C))',Nr,Nr);
 
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invC*B2);
  fp = f1+B1T*(invC*f2);

  % grounding the solution
  S(indc,:) = sparse(1,Np);
  S(:,indc) = sparse(Np,1);
  S(indc,indc) = 1;
  fp(indc)     = 0;

  %
  % partion diagonal and upper diagonal blocks
  %
  diag_block_S       = cell(Nt, 1);
  below_diag_block_S= cell(Nt-1, 1);
  nAi=ncellphi;
  for i=1:Nt-1
    diag_block_S{i}      =S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
    below_diag_block_S{i}=S(i    *nAi+1 : (i+1)*nAi , (i-1)*nAi+1 : i*nAi);
  end
  i=Nt;
  diag_block_S{i}=S((i-1)*nAi+1:i*nAi,(i-1)*nAi+1:i*nAi);

  %				
  % create define inverse on diagonal block
  % relax diagonal by a factor 1/omega 0<omega<1
  %
  diag_block_invS(1:Nt,1)=sparse_inverse;
  for i=1:Nt
    diag_block_invS(i).init((1/ctrl_inner.omega)*diag_block_S{i},ctrl_inner);
  end

  % set initial solution
  x0=zeros(Np,1);
  info_S=info_solver;

  approach="stationary_block_diag"
  %approach="stationary_block_triangular"
  if ( strcmp(approach, "stationary_block_diag")) 
    [dp,info_S.iter,info_S.res,info_S.resvec] = ...
    stationary_iterative_methods(S,...
				 fp,x0,ctrl_outer.tolerance,ctrl_outer.itermax,...
  				 @(x) apply_block_triangular_inverse(...
								      diag_block_invS,below_diag_block_S,'L',x));
    info_S.print();
  elseif ( strcmp(approach, "stationary_block_triangular") )
    [dp,info_S.iter,info_S.res,info_S.resvec] = ...
    stationary_iterative_methods(S,fp,x0,ctrl_outer.tolerance,ctrl_outer.itermax,...
   				 @(x) apply_block_diag_inverse(diag_block_invS,x));
    info_S.print();
    
  end 
 
  % get solution of system JF d  = F  
  dr = invC*(B2*dp-f2);

  %J=[A, B1T;  B2, -C];
  %norm(J*[dp;dr]-[f1;f2])

  ds = JF.ss\(h-JF.sr*dr);
  d = [dp; dr; ds];
    
  normd = norm(d);

elseif sol==9
  % use P=( A BT ) as preconditoner
  %       ( 0 -C )
  % within flexible gmres
  
  
  % grounding
  jacobian(indc,:) = sparse(1,Np+Nr); jacobian(:,indc) = jacobian(Np+Nr,1); jacobian(indc,indc) = 1;
  rhs(indc)=0;

  
  % partion diagonal block of A
  diag_block_A       = cell(Nt, 1);
  nAi=ncellphi;
  for i=1:Nt
    diag_block_A{i}      =A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
  end


 % init A^{-1}=block[A_1^{-1}, A_^{-1}]
  disp('INIT block A^{-1}')
  relax_A=1e-5;
  diag_block_invA(1:Nt,1)=sparse_inverse;
  for i=1:Nt
    diag_block_invA(i).init(diag_block_A{i}+relax_A*speye(nAi),ctrl_inner);
  end

  % init C^{-1}
  fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  invC = sparse(1:Nr,1:Nr,(1.0./spdiags(C))',Nr,Nr);

 % assembly schur AC S=A+B1
 % reduced to system only in phi
  S = A+B1T*(invC*B2);
  fp = f1+B1T*(invC*f2);
  
				% grounding the solution
  S(indc,:) = sparse(1,Np);
  S(:,indc) = sparse(Np,1);
  S(indc,indc) = 1;

  inv_S=sparse_inverse;
  inv_S.init(S,ctrl_inner);

  % solve with fgmres
  stationary=1
  info_S=info_solver;
  if ( stationary)
    [d,info_S.flag,info_S.res,info_S.iter,info_S.resvec] = gmres(jacobian,rhs,20,...
				 ctrl_outer.tolerance,...
				 ctrl_outer.itermax,...
				 @(x) upper_triang_prec(x, ...
							@(y) inv_S.apply(y),...
							B1T,...
							@(z) invC*(-z)));
    	%@(y) apply_block_diag_inverse(diag_block_invA,y) ,...% A^{-1}
    
    info_S.print();
    flag=info_S.flag
    relres=info_S.res
  else
    [d,info_S.flag,info_S.iter,info_S.res,info_S.resvec] = fgmres(jacobian,b,ctrl_outer.tolerance,...
		      'max_iters',ctrl_outer.itermax,...
                      'restart',20,...
		      'verb',0,...
		      'P',@(x,tol) triang_prec(x, ...
					       @(y) apply_block_diag_inverse(diag_block_invA,y) ,...% A^{-1}
					       B1T,...
					       @(z) invC*(-z)));
    info_S.print();
    flag = 0;
    relres = 0;

  end

  res=norm(jacobian*d-rhs)/norm(rhs)

  
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  
  
  sol_stat = struct('flag',flag,'relres',relres,'iter',iter);

elseif sol==10
  verbose=0;
  % Krylov method for full system
  % + 
  % use P=( SAC BT ) as preconditoner 
  %       ( 0   -C )

  % cpu time preprocessing main solver
  tic;
  
  % grounding
  jacobian(indc,:) = sparse(1,Np+Nr); jacobian(:,indc) = jacobian(Np+Nr,1); jacobian(indc,indc) = 1;
  rhs(indc)=0;

  A(indc,:)   = sparse(1,Np);
  A(indc,indc)= 1;
  B1T(indc,:) = sparse(1,Nr);
  B2(:,indc)  = sparse(Nr,1);
  


  
 % init C^{-1}
  if ( verbose >= 1)
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  end
  invC = sparse(1:Nr,1:Nr,(1.0./spdiags(C))',Nr,Nr);

  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invC*B2);
  fp = f1+B1T*(invC*f2);

  if ( compute_eigen)
    eigvalues=eig(full(S));
    fprintf('eig(S) min %8.4e max %8.4e  \n', min(eigvalues), max(eigvalues))
  end
    
  
  
  % grounding the solution
  if ( verbose >= 2)
    fprintf('INIT S=A+B1T C^{-1} B2\n')
  end
 
  %S(indc,:) = sparse(1,Np);
  %S(:,indc) = sparse(Np,1);
  %S(indc,indc) = 1;
  preprocess_cpu=toc;

  approach='full';  
  relax4prec=1e-5;
  debug=0;
  if ( verbose >= 2)
    fprintf('APPROACH %s \n',appraoch)
  end
  if (strcmp(approach,'full'))
    % init inverse of full S
    if ( verbose >= 2)
      fprintf('INIT inverse S\n')
    end
    inv_S=sparse_inverse;
    inv_S.init(S+relax4prec*speye(Np,Np),ctrl_inner);
    inv_S.cumulative_iter=0;
    inv_S.cumulative_cpu=0;
				%inv_S.info();
    inverse_cpu=inv_S.init_cpu;
    if ( verbose >= 2)
      fprintf('DONE inverse S\n') 
    end

    inverse_block11 = @(x) inv_S.apply(x);
    inner_nequ=Np;
  elseif(strcmp(approach,'block_triang'))
    %
    % partion diagonal and upper diagonal blocks
    %
    if ( verbose >= 2)
      fprintf('PARTION diag lower part S\n')
    end
    diag_block_S       = cell(Nt, 1);
    below_diag_block_S = cell(Nt-1, 1);

    if (debug)
      explicit_diag_block_S=sparse(Np,Np);
    end
    
    nAi=ncellphi;
    for i=1:Nt-1
      diag_block_S{i}      =S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      below_diag_block_S{i}=S(i    *nAi+1 : (i+1)*nAi , (i-1)*nAi+1 : i*nAi);
      if (debug)
	explicit_diag_block_S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)=...
	S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
	explicit_diag_block_S(i    *nAi+1 : (i+1)*nAi , (i-1)*nAi+1 : i*nAi)=...
	S(i    *nAi+1 : (i+1)*nAi , (i-1)*nAi+1 : i*nAi);
      end
    end
    i=Nt;
    diag_block_S{i}=S((i-1)*nAi+1:i*nAi,(i-1)*nAi+1:i*nAi);
    if (debug)
      explicit_diag_block_S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)=...
      S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
    end

    
    %				
    % create define inverse on diagonal block
    % relax diagonal by a factor 1/omega 0<omega<1
    %
    if ( verbose >= 2)
      fprintf('INIT INVERSE diag_block(S)\n')
    end 
    diag_block_invS(Nt,1)=sparse_inverse;
    for i=1:Nt
      diag_block_invS(i).init(diag_block_S{i}+relax4prec*speye(nAi,nAi),ctrl_inner);
    end
    
    inverse_block11 = @(x) apply_block_triangular_inverse(diag_block_invS,below_diag_block_S,'L',x);
    if (debug)
      inverse_block11 = @(x) explicit_diag_block_S\x;
    end
    inner_nequ=nAi;
  elseif(strcmp(approach,'block_diag'))
    %
    % partion diagonal and upper diagonal blocks
    %
    if ( verbose >= 2)
      fprintf('PARTION diag S\n')
    end
    diag_block_S       = cell(Nt, 1);
    nAi=ncellphi;
    if (debug) 
      explicit_diag_block_S=sparse(Np,Np);
    end
    for i=1:Nt
      diag_block_S{i}      =S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      if (debug) 
	explicit_diag_block_S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)=...
	S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      end
    end

			%				
			% create define inverse on diagonal block
			% relax diagonal by a factor 1/omega 0<omega<1
			%
    if ( verbose >= 2)
      fprintf('INIT INVERSE diag_block(S)\n')
    end 
    diag_block_invS(Nt,1)=sparse_inverse;
    for i=1:Nt
      diag_block_invS(i).init(diag_block_S{i}+relax4prec*speye(nAi,nAi),ctrl_inner);
    end
    inverse_block11  = @(x) apply_block_diag_inverse(diag_block_invS,x);
    inverse_block11  = @(x) explicit_diag_block_S\x;
    inner_nequ=nAi;
    
  elseif (strcmp(approach,'fullA'))
    % init inverse of full S
    if ( verbose >= 2)
      fprintf('INIT inverse S\n')
    end
    inv_A=sparse_inverse;
    inv_A.init(A+relax4prec*speye(Np,Np),ctrl_inner);
    inv_A.cumulative_iter=0;
    inv_A.cumulative_cpu=0;
				%inv_S.info();
    inverse_cpu=inv_A.init_cpu;
    if ( verbose >= 2)
      fprintf('DONE inverse S\n') 
    end

    inverse_block11 = @(x) inv_A.apply(x);
    inner_nequ=Np;
  end

  %
  % select lower or upper triang
  %
  
  upper_or_lower='lower';
  if (strcmp(upper_or_lower,'lower') )
	    prec = @(x,tol) lower_triang_prec(x, ...
					  @(y) inverse_block11(y),...
					  B2,...
					  @(z) invC*(-z));
  elseif (strcmp(upper_or_lower,'upper') )
	   prec = @(x) upper_triang_prec(x, ...
					 @(y) inverse_block11(y),...
					 B1T,...
					 @(z) invC*(-z));
  end
				% solve with fgmres or bicgstab if prec is linear
  outer_timing=tic;
  info_J=info_solver;
  if ( strcmp(ctrl_inner.approach,'\') | (ctrl_inner.tolerance < 1e-12) )
				%    'USING BICGSTAB'
    [d,info_J.flag,info_J.res,info_J.iter,info_J.resvec] = bicgstab(jacobian,rhs,...
				 ctrl_outer.tolerance,...
				 ctrl_outer.itermax,...
				 @(x) prec(x));    
  else
    [d,i,info_J.resvec,info_J.iter] = fgmres(jacobian,rhs,ctrl_outer.tolerance,...
		      'max_iters',ctrl_outer.itermax,...
                      'restart',20,...
		      'verb',0,...
		      'P',@(x,tol) prec(x));
  end

  if(strcmp(approach,'full'))
    inner_iter=inv_S.cumulative_iter;
    inner_cpu=inv_S.cumulative_cpu;
  elseif(strcmp(approach,'block_diag'))
    inner_iter=0;
    for i=1:Nt
      inner_iter = inner_iter + diag_block_invS(i).cumulative_iter;
    end
    inner_cpu=0.0;
  elseif(strcmp(approach,'block_triang'))
    inner_iter=0;
    for i=1:Nt
      inner_iter = inner_iter + diag_block_invS(i).cumulative_iter;
    end
    inner_cpu=0.0;
  end

  
  outer_iter=uint8(info_J.iter);
  outer_cpu=toc(outer_timing);
  
  info_J.resvec=info_J.resvec(1:outer_iter);
  info_J.res=norm(jacobian*d-rhs)/norm(rhs);
  info_J.flag=1;
  if (info_J.res <= ctrl_outer.tolerance)
    info_J.flag=0;
  end
 
%info_J.print();								      diag_block_invS,below_diag_block_S,'L',y),...
  
				%disp('END')
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  relres=norm(jacobian*d-rhs)/norm(rhs);

  
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

elseif sol==11
  verbose=0;
  % Krylov method for full system
  % + 
  % use P=( ~diag(A)  BT ) as preconditoner 
  %       ( B         -C )

  % cpu time preprocessing main solver
  tic;
  
  % grounding
  jacobian(indc,:) = sparse(1,Np+Nr); jacobian(:,indc) = jacobian(Np+Nr,1); jacobian(indc,indc) = 1;
  rhs(indc)=0;

  A(indc,:)   = sparse(1,Np);
  A(indc,indc)= 1;
  B1T(indc,:) = sparse(1,Nr);
  B2(:,indc)  = sparse(Nr,1);
  


  
 % init C^{-1}
  if ( verbose >= 1)
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  end
  invC     = sparse(1:Nr,1:Nr,(1.0./spdiags(C,0))',Nr,Nr);
  invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);


  % assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
  % reduced to system only in rho
  SCA = C+B2*invDiagA*B1T;
  nnz(SCA)
  fr = f2-B2*(invDiagA*f1);


  % preprocess of precondtioner defienition
  preprocess_cpu=toc;

  approach='full';  
  relax4prec=0e-5;
  debug=0;
  if ( verbose >= 2)
    fprintf('APPROACH %s \n',appraoch)
  end
  
  
  % init inverse of full S
  if ( verbose >= 2)
    fprintf('INIT inverse S\n')
  end
  inv_SCA=sparse_inverse;
  inv_SCA.init(SCA+relax4prec*speye(Nr,Nr),ctrl_inner);
  inv_SCA.cumulative_iter=0;
  inv_SCA.cumulative_cpu=0;
  
  inverse_cpu=inv_SCA.init_cpu;
  if ( verbose >= 2)
    fprintf('DONE inverse S\n') 
  end 
  inverse_block22 = @(x) -inv_SCA.apply(x);
  inner_nequ=Nr;

  % inverse A action
  inverseA_approach='diag';
  %inverseA_approach='full';
  %inverseA_appraoch='block';
  relax4prec=1e-12;
  if ( strcmp(inverseA_approach,'diag'))
    invA = @(x) invDiagA*x;
  elseif( strcmp(inverseA_approach,'full'))
    inv_A=sparse_inverse;
    inv_A.init(A+relax4prec*speye(Np,Np),ctrl_inner);
    inv_A.cumulative_iter=0;
    inv_A.cumulative_cpu=0;

    invA = @(x) inv_A.apply(x);
  elseif( strcmp(inverseA_approach,'block'))
    diag_block_A = cell(Nt, 1);
    nAi=ncellphi;
    for i=1:Nt
      diag_block_A{i}      =A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      if (debug) 
	explicit_diag_block_S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)=...
	S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      end
    end
    diag_block_invA(Nt,1)=sparse_inverse;
    for i=1:Nt
      diag_block_invA(i).init(diag_block_A{i}+relax4prec*speye(nAi,nAi),ctrl_inner);
    end
    invA  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end
  % Define action of preconditoner
  prec = @(x,tol) SchurCA_based_preconditioner(x, invA,inverse_block22,B1T,B2);

  % solve with fgmres or bicgstab if prec is linear
  outer_timing=tic;
  info_J=info_solver;
  if ( strcmp(ctrl_inner.approach,'\') | (ctrl_inner.tolerance < 1e-12) )
    disp( 'USING BICGSTAB')
    [d,info_J.flag,info_J.res,info_J.iter,info_J.resvec] = bicgstab(jacobian,rhs,...
				 ctrl_outer.tolerance,...
				 ctrl_outer.itermax,...
				 @(x) prec(x));    
  else
    [d,i,info_J.resvec,info_J.iter] = fgmres(jacobian,rhs,ctrl_outer.tolerance,...
		      'max_iters',ctrl_outer.itermax,...
                      'restart',20,...
		      'verb',0,...
		      'P',@(x,tol) prec(x));
  end

  inner_iter=inv_SCA.cumulative_iter;
  inner_cpu=inv_SCA.cumulative_cpu;
 
  outer_iter=uint64(info_J.iter);
  outer_cpu=toc(outer_timing);
  
  info_J.resvec=info_J.resvec(1:outer_iter);
  info_J.res=norm(jacobian*d-rhs)/norm(rhs);
  info_J.flag=1;
  if (info_J.res <= ctrl_outer.tolerance)
    info_J.flag=0;
  end
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  relres=norm(jacobian*d-rhs)/norm(rhs);
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

elseif sol==12
  verbose=0;
  % Krylov method for full system
  % + 
  % use P=SchurCA approximation as preconditoner 
  %       

  % cpu time preprocessing main solver
  tic;
  
  % grounding
  jacobian(indc,:) = sparse(1,Np+Nr); jacobian(:,indc) = jacobian(Np+Nr,1); jacobian(indc,indc) = 1;
  rhs(indc)=0;

  A(indc,:)   = sparse(1,Np);
  A(indc,indc)= 1;
  B1T(indc,:) = sparse(1,Nr);
  B2(:,indc)  = sparse(Nr,1);
  
  % init inverse A
  relax4prec=1e-12;
  inverseA_approach='full'
  %inverseA_approach='block'
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    inv_A=sparse_inverse;
    inv_A.init(A+relax4prec*speye(Np,Np),ctrl_inner);
    inv_A.cumulative_iter=0;
    inv_A.cumulative_cpu=0;

    % define function
    invA = @(x) inv_A.apply(x);
  elseif( strcmp(inverseA_approach,'block'))
    % partion matrix 
    diag_block_A = cell(Nt, 1);
    nAi=ncellphi;
    for i=1:Nt
      diag_block_A{i}      =A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
    end
    % use block inverse
    diag_block_invA(Nt,1)=sparse_inverse;
    for i=1:Nt
      diag_block_invA(i).init(diag_block_A{i}+relax4prec*speye(nAi,nAi),ctrl_inner);
    end

    % define function
    invA  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end

  
  % assembly, implicetely, SAC=-(C+B2 * (A)^{-1} B1T)
  SCA = @(x) -(C*x + B2*inv_A.apply(B1T*x))

  % assembly inverse of  SAC using pcg
  ctrl_fmgres=ctrl_solver;
  ctrl_fmgres.init('krylov',1e-5,1000,1.0,1);
  inv_SCA = @(x) fgmres(@(y,tol) SCA(y) ,x,ctrl_fmgres.tolerance,...
		      'max_iters',ctrl_fmgres.itermax,...
                      'restart',20,...
		      'verb',0)
  inner_nequ=Nr;

  % Define action of preconditoner
  prec = @(x,tol) SchurCA_based_preconditioner(x, invA,inv_SCA,B1T,B2);

  
  %preprocess finished
  preprocess_cpu=toc;



  % solve with fgmres or bicgstab if prec is linear
  outer_timing=tic;
  info_J=info_solver;
  [d,i,info_J.resvec,info_J.iter] = fgmres(jacobian,rhs,ctrl_outer.tolerance,...
		      'max_iters',ctrl_outer.itermax,...
                      'restart',20,...
		      'verb',0,...
		      'P',@(x,tol) prec(x));

  
  inner_iter=inv_SCA.cumulative_iter;
  inner_cpu=inv_SCA.cumulative_cpu;
 
  outer_iter=uint64(info_J.iter);
  outer_cpu=toc(outer_timing);
  
  info_J.resvec=info_J.resvec(1:outer_iter);
  info_J.res=norm(jacobian*d-rhs)/norm(rhs);
  info_J.flag=1;
  if (info_J.res <= ctrl_outer.tolerance)
    info_J.flag=0;
  end
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  relres=norm(jacobian*d-rhs)/norm(rhs);
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  
end

realres=compute_linear_system_residuum(JF,F,d);
sol_stat = struct('flag',flag,...
		  'relres',realres,...
		  'ressol',res,...
		  'iter',iter,...
		  'rhsnorm',norm([F.p;F.r;F.s]),...
		  'inner_nequ', inner_nequ,...
		  'inner_iter', inner_iter,...
		  'inner_cpu',  inner_cpu,...
		  'outer_iter', outer_iter,...
		  'outer_cpu',  outer_cpu );
