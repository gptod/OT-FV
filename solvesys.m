function [d, sol_stat, normd] = solvesys(JF,F,controls,logID)

% solve the linear system d = -JF\F, reducing the system by eliminating the
% additional variable s and imposing the condition d(indc) = 0 to fix the 
% constant of the potential dp


indc=controls.indc;
sol=controls.sol;
ctrl_inner11=controls.ctrl_inner11;
ctrl_inner22=controls.ctrl_inner22;
ctrl_outer=controls.ctrl_outer;
compute_eigen=controls.compute_eigen;
  
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

    inner_iter=0;
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
  invS.init(S,ctrl_inner11);
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
 
  % grounding
  A(indc,:)= sparse(1,Np); %A(:,indc)= sparse(Np,1);
  A(indc,indc)=1;
  B1T(indc,:)= sparse(1,Nr);
  rhs(indc)=0;

  
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invdiagC*B2);
  fp = rhs(1:Np)+B1T*(invdiagC*rhs(1+Np:Np+Nr));

  if ( compute_eigen )
    disp("START EIG")
    eigenvalues=eig(full(S));
    disp("END EIG")
    for i=1:4
      fprintf('(%9.4e <= eigvalues <= %18.12e) \n',eigenvalues(i),eigenvalues(Np-i+1))
    end
  end
    
  % solve

  % init solver
  invS=sparse_inverse;
  invS.init(S,ctrl_inner11);
				% apply inverse
  dp=invS.apply(fp);
  inner_nequ=invS.nequ;
  
  % copy info
  flag=invS.info_inverse.flag;
  relres=invS.info_inverse.res;
  inner_iter=invS.info_inverse.iter;
  outer_iter=0;
  outer_cpu=0.0;
  invS.kill();
  
  % get solution of system JF d  = F  
  dr = invdiagC*(B2*dp-rhs(1+Np:Np+Nr));

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

  % grounding
  A(indc,:)= sparse(1,Np);
  A(indc,indc)=1;
  B1T(indc,:)= sparse(1,Nr);
  rhs(indc)=0;

  
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invC*B2);
  fp = rhs(1:Np)+B1T*(invC*rhs(1+Np:Np+Nr));

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
    diag_block_invS(i).init((1/ctrl_inner11.omega)*diag_block_S{i},ctrl_inner11);
  end

  % set initial solution
  x0=zeros(Np,1);
  info_S=info_solver;

  approach=controls.extra_info;
			       %approach="stationary_block_triangular"
  if ( strcmp(approach, 'block_diag')) 
    prec = @(x) apply_block_triangular_inverse(diag_block_invS,below_diag_block_S,'L',x);
  elseif ( strcmp(approach, 'block_triangular') )
    prec = @(x) apply_block_diag_inverse(diag_block_invS,x);
  end

  % solve 
  tic;
  [dp,info_S] = apply_iterative_solver(@(y) S*y, fp, ctrl_outer, prec,x0);
  outer_cpu = toc;
  
  
  % get solution of system JF d  = F  
  dr = invC*(B2*dp-rhs(Np+1:Np+Nr));
  ds = JF.ss\(h-JF.sr*dr);
  d = [dp; dr; ds];

  % store info
  inner_nequ=diag_block_invS(1).nequ;
  inner_iter=0;
  for i  = 1:Nt
    inner_iter=inner_iter+diag_block_invS(i).cumulative_iter;
  end
  outer_iter=info_S.iter

    
  normd = norm(d);

elseif sol==9
  % use P=( A BT ) as preconditoner
  %       ( 0 -C )
  % within flexible gmres
  
  
  % grounding
  %jacobian(indc,:) = sparse(1,Np+Nr); jacobian(:,indc) = jacobian(Np+Nr,1); jacobian(indc,indc) = 1;
  A(indc,:)= sparse(1,Np);
  %A(:,indc)= sparse(Np,1);
  A(indc,indc)=1;
  B1T(indc,:)= sparse(1,Nr);
  rhs(indc)=0;

  
  % partion diagonal block of A
  diag_block_A       = cell(Nt, 1);
  nAi=ncellphi;
  for i=1:Nt
    diag_block_A{i}      =A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
  end


 % init A^{-1}=block[A_1^{-1}, A_^{-1}]
  relax_A=1e-3;
  use_block_invA=1;
  if ( use_block_invA ) 
    diag_block_invA(Nt,1)=sparse_inverse;
    for i=1:Nt
      diag_block_invA(i).init(diag_block_A{i}+relax_A*speye(nAi),ctrl_inner11);
      label=strcat('A_',num2str(i));
      diag_block_invA(i).info_inverse.label=label;
    end
  else
    invA=sparse_inverse;
    invA.init(A+relax_A*speye(Np),ctrl_inner11);
  end
  
  % init C^{-1}
  fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  invC = sparse(1:Nr,1:Nr,(1.0./spdiags(C))',Nr,Nr);

  
  % solve 
  stationary=1;
  info_J=info_solver;
  tic;
  if ( stationary)
    if (use_block_invA)
      [d,info_J.flag,info_J.res,...
       info_J.iter,info_J.resvec] = ...
      bicgstab(jacobian,rhs,...
	       ctrl_outer.tolerance,...
	       ctrl_outer.itermax,...
	       @(x) upper_triang_prec(x, ...
				      @(y) apply_block_diag_inverse(diag_block_invA,y) ,...% A^{-1}
				      B1T,...
				      @(z) -invC*z ));
				% set info solver quantities
      inner_iter=0;
      for i=1:Nt
	inner_iter=inner_iter+diag_block_invA(i).cumulative_iter;
      end
      inner_nequ=nAi;
      outer_iter = info_J.iter;
    else
      [d,info_J.flag,info_J.res,...
       info_J.iter,info_J.resvec] = ...
      bicgstab(jacobian,rhs,...
	       ctrl_outer.tolerance,...
	       ctrl_outer.itermax,...
	       @(x) upper_triang_prec(x, ...
				      @(y) invA.apply(y),...% A^{-1}
				      B1T,...
				      @(z) -invC*z ));

      inner_iter=invA.cumulative_iter;
      inner_nequ=Np;
      outer_iter = info_J.iter;
    end

    flag=info_J.flag;
    relres=info_J.res;
  else
    [d,info_J.flag,iters,info_J.res,info_J.resvec] = fgmres(jacobian,rhs,ctrl_outer.tolerance,...
		      'max_iters',ctrl_outer.itermax,...
                      'restart',20,...
		      'verb',0,...
		      'P',@(x,tol) triang_prec(x, ...
					       @(y) apply_block_diag_inverse(diag_block_invA,y) ,...% A^{-1}
					       B1T,...
					       @(z) invC*(-z)));
    info_J.iter=iters(1);
    info_J.print();
    flag = 0;
    relres = 0;

  end
  outer_cpu=toc;

  % compute real residuum
  A = sparse(JF.pp); B1T = sparse(JF.pr); B2 = sparse(JF.rp);
  C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
  jacobian2 = [A B1T; B2 -C];
  relres=norm(jacobian2*d-rhs)/norm(rhs);

  % 
  
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  
  
  sol_stat = struct('flag',flag,'relres',relres,'iter',iter);

elseif sol==10
  verbose=0;
  % Krylov method for full system
  % + 
  % use P=( SAC B1T ) as preconditoner 
  %       ( 0    -C )

  % cpu time preprocessing main solver
  tic;
  
  % grounding
  if (indc>0)
    rhs(indc)=0;  
    A(indc,:)   = sparse(1,Np);
    A(indc,indc)= 1;
    B1T(indc,:) = sparse(1,Nr);
    B2(:,indc)  = sparse(Nr,1);
  end


  
   % init diag(C)^{-1}
  invC = sparse_inverse;
  invC.init(C,ctrl_inner22);
  invdiagC = sparse(1:Nr,1:Nr,(invC.inverse_matrix_diagonal)',Nr,Nr);
  if ( verbose >= 1)
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
  end

  
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invdiagC*B2);
  %fp = rhs(1:Np)+B1T*(invC*rhs(Np+1:Np+Nr));

  if ( compute_eigen)
    eigvalues=eig(full(S));
    fprintf('eig(S) min %8.4e max %8.4e  \n', min(eigvalues), max(eigvalues))
  end
  
  % grounding the solution
  if ( verbose >= 2)
    fprintf('INIT S=A+B1T C^{-1} B2\n')
  end
 

  % assembly (approcimate) inverse of SAC=A+B1T diag(C)^{-1} B2 
  preprocess_cpu=toc;
  approach=controls.extra_info;  
  relax4prec=controls.relax4inv11;
  relax4prec=0.0;
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
    inv_S.init(S+relax4prec*speye(Np,Np),ctrl_inner11);
    inv_S.cumulative_iter=0;
    inv_S.cumulative_cpu=0;
    inv_S.dimblock=ncellphi;
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
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
		    ctrl_inner11.tolerance,...
		    ctrl_inner11.itermax,...
		    ctrl_inner11.omega,...
		    ctrl_inner11.verbose,...
		    sprintf('%sA%d%d',ctrl_inner11.label,i,i));
      diag_block_invA(i).name=sprintf('inverse A%d%d',i,i);
      diag_block_invS(i).init(diag_block_S{i}+relax4prec*speye(nAi,nAi),ctrl_loc);
    end
    for i=1:Nt
      disp(diag_block_invA(i).ctrl.label)
    end
      
      
    
    
    inverse_block11 = @(x) apply_block_triangular_inverse(diag_block_invS,below_diag_block_S,'L',x);
    if (debug)
      inverse_block11 = @(x) explicit_diag_block_S\x;
    end
    inner_nequ=nAi;
  elseif(strcmp(approach,'block_diag'))
    nAi=ncellphi;
    inner_nequ=nAi;
    
    %
    % partion diagonal and upper diagonal blocks
    %
    if ( verbose >= 2)
      fprintf('PARTION diag S\n')
    end
    diag_block_S       = cell(Nt, 1);

    
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
    ctrl_loc=ctrl_solver;
    for i=1:Nt
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
		    ctrl_inner11.tolerance,...
		    ctrl_inner11.itermax,...
		    ctrl_inner11.omega,...
		    ctrl_inner11.verbose,...
		    sprintf('%s_SAC%d%d',ctrl_inner11.label,i,i));
      diag_block_invS(i).init(diag_block_S{i}+relax4prec*speye(nAi,nAi),ctrl_loc);
      diag_block_invS(i).name=sprintf('inverse_SAC%d%d',ctrl_inner11.label,i,i);
    end
    if (1) 
      for i=1:Nt
	disp(diag_block_invS(i).ctrl.label);
	Sii=S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
	differe= norm(nonzeros(diag_block_invS(i).matrix)-...
		      nonzeros(Sii));
	fprintf('diff block %1.4e\n', differe);
      end
      inverse_block11  = @(x) apply_block_diag_inverse(diag_block_invS,x);
    end
    
  end

  %
  % Define action of preconditoner
  % 
  prec = @(x) SchurAC_based_preconditioner(x, inverse_block11, @(y) invC.apply(y) ,B1T,B2,controls.outer_prec,ncellphi);

  %prec = @(x) lower_triang_prec(x,inverse_block11,B2, @(y) invC.apply(y));


  % solve with fgmres or bicgstab if prec is linear
  outer_timing=tic;
  jacobian = [A B1T; B2 -C];
  [d,info_J]=apply_iterative_solver(@(x) jacobian*x, rhs, ctrl_outer, prec );
  outer_cpu=toc(outer_timing);

  if(strcmp(approach,'full'))
    inner_iter=inv_S.cumulative_iter;
  elseif(strcmp(approach,'block_diag'))
    inner_iter=0;
    for i=1:Nt
      inner_iter = inner_iter + diag_block_invS(i).cumulative_iter;
    end
  elseif(strcmp(approach,'block_triang'))
    inner_iter=0;
    for i=1:Nt
      inner_iter = inner_iter + diag_block_invS(i).cumulative_iter;
    end
  end

  
  outer_iter=uint64(info_J.iter);
  
  
  %info_J.resvec=info_J.resvec(1:outer_iter);
  info_J.res=norm(jacobian*d-rhs)/norm(rhs);
  info_J.flag=1;
  if (info_J.res <= ctrl_outer.tolerance)
    info_J.flag=0;
  end
 
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
  verbose=controls.verbose;
  % Krylov method for full system
  % + 
  % use P=( ~A  B1T ) as preconditoner 
  %       ( B2          -C )
  %
  % 
  % based on
  % ~A=diag(A)
  % ~A=approximate inverse of A given by inner_solver
  % S=-C-B2*(diag(A))^{-1} B1T 

  % cpu time preprocessing main solver
  tic;

  % grounding
  if ( indc >0) 
    rhs(indc)=0;
    A(indc,:)   = sparse(1,Np);
    %A(:,indc)   = sparse(Np,1);
    A(indc,indc)= 1;
    B1T(indc,:) = sparse(1,Nr);
    %B2(:,indc) = sparse(Nr,1);
  end
  

  if ( verbose >= 2)
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
    fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)))
  end
  
  % diag(A)^{-1}
  relax4inv11=1e-12;
  inv_A=sparse_inverse;
  inv_A.init(A+relax4inv11*speye(Np,Np),ctrl_inner11);
  inv_A.info_inverse.label='schur_ca';
  inv_A.cumulative_iter=0;
  inv_A.cumulative_cpu=0;

  % define inverse of (~A)^{-1}
  relax4prec=controls.relax4inv11;
  
  % inverse A action
  inverseA_approach=controls.extra_info;
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    invA=sparse_inverse;
    invA.init(A+relax4prec*speye(Np,Np),ctrl_inner11);
    invA.cumulative_iter=0;
    invA.cumulative_cpu=0;

    % define function
    inv_A = @(x) invA.apply(x);
  elseif( strcmp(inverseA_approach,'block'))
    % partion matrix 
    diag_block_A = cell(Nt, 1);
    nAi=ncellphi;
    for i=1:Nt
      diag_block_A{i}      =A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
    end
    % use block inverse
    diag_block_invA(Nt,1)=sparse_inverse;

    ctrl_loc=ctrl_solver;
    for i=1:Nt
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
		    ctrl_inner11.tolerance,...
		    ctrl_inner11.itermax,...
		    ctrl_inner11.omega,...
		    ctrl_inner11.verbose,...
		    sprintf('%sA%d',ctrl_inner11.label,i));
      diag_block_invA(i).name=sprintf('inverse A%d',i);
      diag_block_invA(i).init(diag_block_A{i}+relax4prec*speye(nAi,nAi),ctrl_loc);     
    end

    % define function
    inv_A  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end

  % define explicitely inverse of diag(A)
  invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);

  
  % assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
  % reduced to system only in rho
  SCA = (C+B2*invDiagA*B1T);
  fr = rhs(1+Np:Np+Nr)-B2*(invDiagA*rhs(1:Np));


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
  inv_SCA.init(SCA+controls.relax4inv22*speye(Nr,Nr),ctrl_inner22);
  inv_SCA.info_inverse.label='schur_ca';
  inv_SCA.cumulative_iter=0;
  inv_SCA.cumulative_cpu=0;
  
  inverse_cpu=inv_SCA.init_cpu;
  if ( verbose >= 2)
    fprintf('DONE inverse S\n') 
  end 
  inverse_block22 = @(x) -inv_SCA.apply(x);

  %SCA = @(x) (C*x + B2*(invA(B1T*x)));
  %assembly inverse of  SAC iterative solver with no precodnitioner
  %inverse_block22 = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);
  inner_nequ=Nr;

  
  % Define action of preconditoner
  prec = @(x) SchurCA_based_preconditioner(x, inv_A,inverse_block22,B1T,B2,controls.outer_prec,ncellphi);

  % solve
  jacobian = [A B1T; B2 -C];

  
  
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) mxv_jacobian(A,B1T,B2,C,ncellphi,x), ...
				    rhs, ctrl_outer, prec,[],controls.left_right );
  outer_cpu=toc(outer_timing);

  % get info
  inner_iter=inv_SCA.cumulative_iter;
  outer_iter=uint32(info_J.iter);
  
  
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  A = sparse(JF.pp); B1T = sparse(JF.pr); B2 = sparse(JF.rp);
  C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
  f1 = f;
  f2 = g-JF.rs*(JF.ss\h);
  % swap C sign for having standard saddle point notation 
  C = -C;


  % assembly full system
  jacobian = [A B1T; B2 -C];

  % assembly rhs
  rhs=[f1;f2];

  

  relres=norm(jacobian*d-rhs)/norm(rhs);
 
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

elseif sol==12
  verbose=controls.verbose;
  % Krylov method for full system
  % + 
  % use P=( ~A  B1T ) as preconditoner 
  %       ( B2  -C )
  %
  % 
  % based on
  % ~A=diag(A)
  % ~A=approximate inverse of A given by inner_solver
  % S=-C-B2*(diag(A))^{-1} B1T 

  % cpu time preprocessing main solver
  tic;
  
				% grounding
  if (indc>0) 
    rhs(indc)=0;
    A(indc,:)   = sparse(1,Np);
    %A(:,indc)   = sparse(Np,1);
    A(indc,indc)= 1;
    B1T(indc,:) = sparse(1,Nr);
    %B2(:,indc) = sparse(Nr,1);
  end
  

  if ( verbose >= 2)
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
    fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)))
  end

  % define inverse of (~A)^{-1}
  relax4prec=1e-8
  
  % inverse A action  
  if ( strcmp(ctrl.extra_info,'full'))
    % set inverse
    invA=sparse_inverse;
    invA.init(A+relax4prec*speye(Np,Np),ctrl_inner11);
    invA.cumulative_iter=0;
    invA.cumulative_cpu=0;

    % define function
    inv_A = @(x) invA.apply(x);
  elseif( strcmp(inverseA_approach,'block'))
    % partion matrix 
    diag_block_A = cell(Nt, 1);
    nAi=ncellphi;
    for i=1:Nt
      diag_block_A{i}      =A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
    end
    % use block inverse
    diag_block_invA(Nt,1)=sparse_inverse;

    ctrl_loc=ctrl_solver;
    for i=1:Nt
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
		    ctrl_inner11.tolerance,...
		    ctrl_inner11.itermax,...
		    ctrl_inner11.omega,...
		    ctrl_inner11.verbose,...
		    sprintf('%s A%d',ctrl_inner11.label,i));
      diag_block_invA(i).name=sprintf('inverse A%d',i);

      
      diag_block_invA(i).init(diag_block_A{i}+relax4prec*speye(nAi,nAi),ctrl_loc);     
    end

    % define function
    inv_A  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end
  
 
  
  % preprocess of precondtioner defienition
  preprocess_cpu=toc;

  % 
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
  

  SCA = @(x) (C*x + B2*(invA(B1T*x)));
  %assembly inverse of  SAC iterative solver with no precodnitioner
  inverse_block22 = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);
  inner_nequ=Nr;

  
  % Define action of preconditoner
  prec = @(x) SchurCA_based_preconditioner(x, invA,inverse_block22,B1T,B2);

  % solve
  jacobian = [A B1T; B2 -C];
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) jacobian*x, rhs, ctrl_outer, prec );
  outer_cpu=toc(outer_timing);

  % get info
  inner_iter=0;%inv_SCA.cumulative_iter;
  outer_iter=uint32(info_J.iter);
  
  
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  A = sparse(JF.pp); B1T = sparse(JF.pr); B2 = sparse(JF.rp);
  C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
  f1 = f;
  f2 = g-JF.rs*(JF.ss\h);
  % swap C sign for having standard saddle point notation 
  C = -C;


  % assembly full system
  jacobian = [A B1T; B2 -C];

  % assembly rhs
  rhs=[f1;f2];

  

  relres=norm(jacobian*d-rhs)/norm(rhs);
 
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);
  
elseif sol==13
  verbose=0;
  % Krylov method for full system
  % + 
  % use P=SchurCA approximation as preconditoner 
  %       

  % cpu time preprocessing main solver
  tic;
  
  % grounding
  rhs(indc)=0;
  A(indc,:)   = sparse(1,Np);
  A(indc,indc)= 1;
  B1T(indc,:) = sparse(1,Nr);

  jacobian = [A B1T; B2 -C];

  % init inverse A
  relax4prec=0e-6;
  inverseA_approach='full';
  %inverseA_approach='block';
  %inverseA_approach='diag';
  if ( strcmp(inverseA_approach,'diag'))
    invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);
    inv_A = @(x) invDiagA*x;
  elseif ( strcmp(inverseA_approach,'full'))
    % set inverse
    invA=sparse_inverse;
    invA.init(A+relax4prec*speye(Np,Np),ctrl_inner11);
    invA.cumulative_iter=0;
    invA.cumulative_cpu=0;

    % define function
    inv_A = @(x) invA.apply(x);
  elseif( strcmp(inverseA_approach,'block'))
    % partion matrix 
    diag_block_A = cell(Nt, 1);
    nAi=ncellphi;
    for i=1:Nt
      diag_block_A{i}      =A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
    end
    % use block inverse
    diag_block_invA(Nt,1)=sparse_inverse;

    ctrl_loc=ctrl_solver;
    disp(' ')
    for i=1:Nt
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
		    ctrl_inner11.tolerance,...
		    ctrl_inner11.itermax,...
		    ctrl_inner11.omega,...
		    ctrl_inner11.verbose,...
		    sprintf('%s A%d',ctrl_inner11.label,i));
      diag_block_invA(i).name=sprintf('inverse A%d',i);

      
      %disp(i)
      %disp(diag_block_invA(i).name)
      diag_block_invA(i).init(diag_block_A{i}+relax4prec*speye(nAi,nAi),ctrl_loc);
      % if (i>1)
      % 	disp(strcat(num2str(i),diag_block_invA(i).ctrl.label,...
      % 		    'before',diag_block_invA(i-1).name,diag_block_invA(i-1).ctrl.label))
      % end
      
    end

    % for i=1:Nt
    %   fprintf( '%d %s %s\n', i,diag_block_invA(i).name,diag_block_invA(i).ctrl.label)	    
    % end


    % define function
    inv_A  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end

  %f_rhs=rhs(1:Np);
  %test=inv_A(f_rhs);
  %resf=norm(A*test-f_rhs)/norm(f_rhs)

  % schur CA
  
  % %assembly, implicetely, SAC=-(C+B2 * (A)^{-1} B1T)
  % SCA = @(x) (C*x + B2*(inv_A(B1T*x)));
  % %assembly inverse of  SAC iterative solver with no precodnitioner
  % invS = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);

  
  invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);
  SCA = (C+B2*invDiagA*B1T);

  relax4prec=0.0
  inv_SCA=sparse_inverse;
  inv_SCA.init(SCA+relax4prec*speye(Nr,Nr),ctrl_inner22);
  inv_SCA.info_inverse.label='schur_ca';
  inv_SCA.cumulative_iter=0;
  inv_SCA.cumulative_cpu=0;
  
  inverse_cpu=inv_SCA.init_cpu;
  if ( verbose >= 2)
    fprintf('DONE inverse S\n') 
  end 
  inverse_block22 = @(x) -inv_SCA.apply(x);
  inner_nequ=Nr;

  

  % Define action of preconditoner
  prec = @(x,tol) SchurCA_based_preconditioner(x, inv_A,inverse_block22,B1T,B2);

  %preprocess finished
  preprocess_cpu=toc;


  % solve 
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) jacobian*x, rhs, ctrl_outer, prec );
  outer_cpu=toc(outer_timing);
  
  flag=info_J.flag;
  relres=info_J.res;
  outer_iter=info_J.iter;


  inner_iter=0;
  inner_cpu=0;
  inner_nequ=0;
  
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
		  'outer_iter', outer_iter,...
		  'outer_cpu',  outer_cpu );
