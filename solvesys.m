function [d, sol_stat, normd] = solvesys(JF,F,controls,logID)

% solve the linear system d = -JF\F, reducing the system by eliminating the
% additional variable s and imposing the condition d(indc) = 0 to fix the 
% constant of the potential dp


sol=controls.sol;
indc=controls.indc;

%ctrl_inner11=controls.ctrl_inner11;
%ctrl_inner22=controls.ctrl_inner22;
%ctrl_outer=controls.ctrl_outer;
%compute_eigen=controls.compute_eigen;
  
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

% assembly rhs
rhs=[f1;f2];

kernel=ones(ncellphi,1);
kernel=kernel/norm(kernel);

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
  %invS.kill();

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
  %invS.kill();
  
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
  % copy controls
  indc=controls.indc;
  sol=controls.sol;
  ctrl_inner11=controls.ctrl_inner11;
  ctrl_inner22=controls.ctrl_inner22;
  ctrl_outer=controls.ctrl_outer;
  compute_eigen=controls.compute_eigen;


  verbose=0;
  % Krylov method for full system
  % + 
  % use P=( SAC B1T ) as preconditoner 
  %       ( 0    -C )

  % cpu time preprocessing main solver
  tic;

  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls,controls.manipulation_approach);
  if (swap_sign)
    vectors_x=-vectors_x;
    vectors_y=-vectors_y;
    alphas=-alphas;
  end
  for i=1:N
    fprintf('(%1.4e <= p1 <= %1.4e) \n',min(vectors_x(:,1+(i-1)*2)),max(vectors_x(:,1+(i-1)*2)))
    fprintf('(%1.4e <= p2 <= %1.4e) \n',min(vectors_x(:,2+(i-1)*2)),max(vectors_x(:,2+(i-1)*2)))
    fprintf('A(%d,%d)=%1.4e \n', indeces_global(i),indeces_global(i),...
   	    find(A(indeces_global(i+1),indeces_global(i+1))));
  end
  
  if (controls.manipulate > 0 )
    [A,B1T,rhs]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas);
  end
  
 % grounding
  if (indc>0)
    % find nodes at time time step diag(A^i) is max
    disp('GROUNDING')
    indeces=set_grounding_node(A,ncellphi);
    inode=indeces(1);
    [A,B1T,rhs]=grounding(A,B1T,rhs,inode,0);
  elseif(indc<0)
    irow=indeces_global(1);
    A(irow,1:ncellphi) = -vectors_x(:,1)';
    B1T(irow,:)=sparse(1,Nr);
    rhs(irow,:)=0.0;
  end


  
   % init diag(C)^{-1}
  invC = sparse_inverse;

  ctrl_loc=ctrl_solver;
  index_agmg=1;
  ctrl_loc.init(ctrl_inner22.approach,...
		ctrl_inner22.tolerance,...
		ctrl_inner22.itermax,...
		ctrl_inner22.omega,...
		ctrl_inner22.verbose,...
		'C',index_agmg);
 
  invC.init(C,ctrl_loc);
  invdiagC = sparse(1:Nr,1:Nr,(invC.inverse_matrix_diagonal)',Nr,Nr);
  if ( verbose >= 1)
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C,0)),max(spdiags(C,0)))
  end

  
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(invdiagC*B2);  
  fp = rhs(1:Np)+B1T*(invC.apply(rhs(Np+1:Np+Nr)));

  
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
    fprintf('APPROACH %s \n',approach)
  end
  if (strcmp(approach,'full'))
    % init inverse of full S
    if ( verbose >= 2)
      fprintf('INIT inverse S\n')
    end
    inv_S=sparse_inverse;

    ctrl_loc=ctrl_solver;
    index_agmg=index_agmg+1;
    ctrl_loc.init(ctrl_inner11.approach,...
		  ctrl_inner11.tolerance,...
		  ctrl_inner11.itermax,...
		  ctrl_inner11.omega,...
		  ctrl_inner11.verbose,...
		  'SAC',index_agmg);
    inv_S.init(S+relax4prec*speye(Np,Np),ctrl_loc);
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
      index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
		    ctrl_inner11.tolerance,...
		    ctrl_inner11.itermax,...
		    ctrl_inner11.omega,...
		    ctrl_inner11.verbose,...
		    sprintf('%sA%d%d',ctrl_inner11.label,i,i),index_agmg);
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
    % create define inverse on diagonal block
    % relax diagonal by a factor 1/omega 0<omega<1
    %
    if ( verbose >= 2)
      fprintf('INIT INVERSE diag_block(S)\n')
    end 
    diag_block_invS(Nt,1)=sparse_inverse;
    ctrl_loc=ctrl_solver;
    for i=1:Nt
      index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
		    ctrl_inner11.tolerance,...
		    ctrl_inner11.itermax,...
		    ctrl_inner11.omega,...
		    ctrl_inner11.verbose,...
		    sprintf('%s_SAC%d%d',ctrl_inner11.label,i,i),index_agmg);
      diag_block_invS(i).init(S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)+...
			      relax4prec*speye(nAi,nAi),ctrl_loc);
      diag_block_invS(i).name=sprintf('inverse_SAC%d%d',ctrl_inner11.label,i,i);
    end

    
    inverse_block11  = @(x) apply_block_diag_inverse(diag_block_invS,x);
  end

  if( 0 )
    sol_phi=inverse_block11(fp);
    sol_rho=invC.apply(B2*sol_phi-rhs(Np+1:Np+Nr));
    fprintf(' res |S x -f_r|/f_r = %1.4e \n',norm(S*sol_phi-fp)/norm(fp));
    
    
    jacobian=[A B1T;B2 -C];
    d=[sol_phi;sol_rho];
    
    
    fprintf(' res modified = %1.4e \n',norm(jacobian*d-rhs)/norm(rhs))
  end
  
  
  %
  % Define action of preconditoner
  % 
  prec = @(x) SchurAC_based_preconditioner(x, inverse_block11, @(y) invC.apply(y) ,...
					   B1T,B2,controls.outer_prec,ncellphi);

  %prec = @(x) lower_triang_prec(x,inverse_block11,B2, @(y) invC.apply(y));


  % solve with fgmres or bicgstab if prec is linear
  outer_timing=tic;
  jacobian = [A B1T; B2 -C];
  if ( verbose >= 2)
    fprintf('START SOLVER \n')
  end
  [d,info_J]=apply_iterative_solver(@(x) jacobian*x, rhs, ctrl_outer, prec );
  if ( verbose >= 2)
    fprintf('END SOLVER \n')
  end

  outer_iter=uint64(info_J.iter);
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
 
  % free memory (mandatory for agmg)
  if (strcmp(approach,'full'))
    inv_S.kill();
  elseif(strcmp(approach,'block_diag'))
    for i=1:Nt
      diag_block_invS(i).kill()
    end
  elseif(strcmp(approach,'block_triang'))
    for i=1:Nt
      diag_block_invS(i).kill()
    end
  end
    
 
  if (1)
    sizecell=-spdiags(JF.rs(1:ncellrho,1:ncellrho),0);
      
    % call subroutine to check if the solution satisfies
    % p^t x + q^t y=alpha
    test_vectors(d,vectors_x,vectors_y,alphas,N,Np,Nr)
    for i=1:N
      fprintf(' %d imbalance=%1.2e \n', i, sizecell'*d(Np+1+(i-1)*ncellrho:Np+i*ncellrho))
    end
  end

  
  
  %info_J.resvec=info_J.resvec(1:outer_iter);
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

elseif sol==11
  % copy controls
  indc=controls.indc;
  sol=controls.sol;
  ctrl_inner11=controls.ctrl_inner11;
  ctrl_inner22=controls.ctrl_inner22;
  ctrl_outer=controls.ctrl_outer;
  compute_eigen=controls.compute_eigen;

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


  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls,controls.manipulation_approach);
  if (~swap_sign)
    vectors_x=-vectors_x;
    vectors_y=-vectors_y;
    alphas=-alphas;
  end

  if ( indc >0)
    % find nodes at time time step diag(A_ii) is max
    indeces=set_grounding_node(A,ncellphi);

    % ground node of A_11
    inode=indeces_globals(1);
    [A,B1T,rhs]=grounding(A,B1T,rhs,inode,0);
  end

  if (controls.manipulate)
    disp('MAKING A invertible')
    [A,B1T,rhs]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas);
  end


  if (controls.diagonal_scaling)
    diagA_scaling = sparse(1:Np,1:Np,(1.0./sqrt(spdiags(A,0)))',Np,Np);
    diagC_scaling=sparse(1:Nr,1:Nr,(1.0./sqrt(spdiags(C,0)))',Nr,Nr);

    
    %fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C,0)),max(spdiags(C,0)))
    %fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)))

   

    
    A=  diagA_scaling*A*  diagA_scaling;
    B1T=diagA_scaling*B1T*diagC_scaling;
    B2= diagC_scaling*B2* diagA_scaling;
    C=  diagC_scaling*C*  diagC_scaling;

    %fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C,0)),max(spdiags(C,0)))
    %fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)))

    rhs(1:Np)      =diagA_scaling*rhs(1:Np);
    rhs(1+Np:Np+Nr)=diagC_scaling*rhs(1+Np:Np+Nr);   
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
    % % partion matrix 
    % nAi=ncellphi;
   
    % % use block inverse
    % diag_block_invA(Nt,1)=sparse_inverse;
    % ctrl_loc=ctrl_solver;


    % [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
   
    % for i=1:Nt
    %   ctrl_loc=ctrl_solver;
    %   ctrl_loc.init(ctrl_inner11.approach,...
    % 		    ctrl_inner11.tolerance,...
    % 		    ctrl_inner11.itermax,...
    % 		    ctrl_inner11.omega,...
    % 		    ctrl_inner11.verbose,...
    % 		    sprintf('%sA%d',ctrl_inner11.label,i));
    %   diag_block_invA(i).name=sprintf('inverse A%d',i);
    %   matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + relax4prec*speye(nAi,nAi) ;
    %   inode=indeces_local(i)
    %   matrixAi(inode,:)=sparse(ncellphi,1);
    %   diag_block_invA(i).init(matrixAi,ctrl_loc);     
    % end

    % % define function
    % inv_A  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end

  % define explicitely inverse of diag(A)
  invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);

  inv_A  = @(x) invDiagA*x;
  
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

  
  %fprintf(' %1.4e <=diag(S)< %1.4e\n ', min(inv_SCA.inverse_matrix_diagonal),max(inv_SCA.inverse_matrix_diagonal));

  %SCA = @(x) (C*x + B2*(invA(B1T*x)));
  %assembly inverse of  SAC iterative solver with no precodnitioner
  %inverse_block22 = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);
  inner_nequ=Nr;

  
				% Define action of preconditoner
  
  prec = @(x) SchurCA_based_preconditioner(x, inv_A,...
					   inverse_block22,...
					   @(y) B1T*y,...
					   @(z) B2*z,...
					   Np,Nr,...
					   controls.outer_prec,ncellphi);
  
  % solve
  jacobian = [A B1T; B2 -C];

  
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) mxv_jacobian(x,A,B1T,B2,C,ncellphi,kernel,rhs), ...
				    rhs, ctrl_outer, prec,[],controls.left_right );
  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);  
  end
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
  % copy controls
  indc=controls.indc;
  sol=controls.sol;
  ctrl_inner11=controls.ctrl_inner11;
  ctrl_inner22=controls.ctrl_inner22;
  ctrl_outer=controls.ctrl_outer;
  compute_eigen=controls.compute_eigen;
  
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

  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls,controls.manipulation_approach);
  if (swap_sign)
    vectors_x=-vectors_x;
    vectors_y=-vectors_y;
    alphas=-alphas;
  end
  %for i=1:N
  %  fprintf('(%1.4e <= p1 <= %1.4e) \n',min(vectors_x(:,1+(i-1)*2)),max(vectors_x(:,1+(i-1)*2)))
  %  fprintf('(%1.4e <= p2 <= %1.4e) \n',min(vectors_x(:,2+(i-1)*2)),max(vectors_x(:,2+(i-1)*2)))
  %  fprintf('A(irow,irow)=%1.4e \n', find(A(indeces_global,indeces_global)));
  %end
  
  if (controls.manipulate > 0)
    [A,B1T,rhs]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas);
  end
  				% grounding
  if (indc>0)
    % find nodes at time time step diag(A^i) is max
    disp('GROUNDING')
    inode=indeces_global(1);
    [A,B1T,rhs]=grounding(A,B1T,rhs,inode,0);
  elseif(indc<0)
    irow=indeces_global(1);
    A(irow,1:ncellphi) = -vectors_x(:,1)';
    B1T(irow,:)=sparse(1,Nr);
    rhs(irow,:)=0.0;
  end


  if (controls.diagonal_scaling)
    diagA_scaling = sparse(1:Np,1:Np,(1.0./sqrt(spdiags(A,0)))',Np,Np);
    diagC_scaling=sparse(1:Nr,1:Nr,(1.0./sqrt(spdiags(C,0)))',Nr,Nr);
    
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
    fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)));
    
    A=  diagA_scaling*A*  diagA_scaling;
    B1T=diagA_scaling*B1T*diagC_scaling;
    B2= diagC_scaling*B2* diagA_scaling;
    C=  diagC_scaling*C*  diagC_scaling;

    %fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
    %fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)))

    rhs(1:Np)      =diagA_scaling*rhs(1:Np);
    rhs(1+Np:Np+Nr)=diagC_scaling*rhs(1+Np:Np+Nr);   
  end

  

  % define inverse of (~A)^{-1}
  relax4prec=controls.relax4inv11;

  index_agmg=0;
  
  % inverse A action
  approach_inverse_A=controls.approach_inverse_A;
  if (strcmp(approach_inverse_A,'full'))
    % we need a new seed for agmg, in case is used
    index_agmg=index_agmg+1;
    ctrl_loc=ctrl_solver;
    ctrl_loc.init(controls.ctrl_inner11.approach,...
		  controls.ctrl_inner11.tolerance,...
		  controls.ctrl_inner11.itermax,...
		  controls.ctrl_inner11.omega,...
		  controls.ctrl_inner11.verbose,...
		  '~A',index_agmg);

    
    inverseA=sparse_inverse;
    inverseA.init(A+relax4prec*speye(Np,Np),ctrl_loc);
    inverseA.info_inverse.label='A';
    inverseA.cumulative_iter=0;
    inverseA.cumulative_cpu=0;

    inv_A  = @(x) inverseA.apply(x);
    
  elseif(strcmp(approach_inverse_A,'block'))
    % define array of sparse inverse
    diag_block_invA(Nt,1)=sparse_inverse;

				% create blocks
    nAi=ncellphi;
    ctrl_loc=ctrl_solver;
    for i=1:Nt
      % we need a new seed for agmg, in case is used
      index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
      % passing index_agmg we will use agmg-i 
      ctrl_loc.init(controls.ctrl_inner11.approach,...
		    controls.ctrl_inner11.tolerance,...
		    controls.ctrl_inner11.itermax,...
		    controls.ctrl_inner11.omega,...
		    controls.ctrl_inner11.verbose,...
		    sprintf('%s A%d',ctrl_inner11.label,i),...
		    index_agmg);
      diag_block_invA(i).name=sprintf('inverse A%d',i);
      
      % get block add relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)+relax4prec*speye(nAi,nAi);
      
				% define inverse
      diag_block_invA(i).init( matrixAi,ctrl_loc);     

    end
% define function A shifter( inverse A grounded ( rhs set to zero in some nodes )
    inv_A  = @(x) apply_inverseAtilde(diag_block_invA,A,indeces_global,indeces_local,x);

  end

  
  % preprocess of precondtioner defienition
  preprocess_cpu=toc;

  % compute eigenvalue of schur complement
  if (controls.compute_eigen)
    
    ctrl_loc=ctrl_solver;
    ctrl_loc.init('direct',1e-14,0,0.0,0,'accurate A inverse');
    real_inverseA=sparse_inverse;
    real_inverseA.init(A+relax4prec*speye(Np,Np),ctrl_inner11);
    dense_S=form_minus_SchurCA(@(x) real_inverseA.apply(x),B1T,B2,C);
    
    fprintf('\n')
    eigenvalues=eig(dense_S);
    real_part=real(eigenvalues);
    imag_part=imag(eigenvalues);
    abs_eig=abs(eigenvalues);
    fprintf('max/min=%1.4e | %1.4e<=real(eig(S))<=%1.4e %1.4e<=img(eig(S))<=%1.4e %1.4e<=abs(eig(S))<=%1.4e\n',max(abs_eig)/min(abs_eig),...
	    min(real_part),max(real_part),...
	    min(imag_part),max(imag_part),...
	    min(abs_eig),max(abs_eig))
    
    figure
    plot(eigenvalues,'o')
  
      
    if (strcmp(approach_inverse_A,'full'))
      inverseA.cumulative_application=0;
    elseif (strcmp(approach_inverse_A,'full'))
      diag_block_invA(1,1).cumulative_application=0;
    end

    %dS=diag(dense_S);
    %fprintf('(%9.4e <= diag(S) <= %18.12e) \n',min(dS),max(dS))

    %sqrt_dS=sparse(1:Nr,1:Nr,(1.0./sqrt(dS))',Nr,Nr);

    %new_dense=sqrt_dS*dense_S*sqrt_dS;
    %new_dense=dense_S;
    %new_dense(abs(new_dense)<1e-10)  = 0;
    
    
    %figure
    %spy(new_dense)
    %imagesc(new_dense)
  end
    

  % 
  approach_schurCA=controls.approach_schurCA;
  if (strcmp(approach_schurCA,'diagA'))
    % define explicitely inverse of diag(A)
    invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);

    % assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
    approx_SCA = (C+B2*invDiagA*B1T);

    % assembly approximate inverse
    inv_SCA=sparse_inverse;


    % we need a new seed for agmg, in case is used
    index_agmg=index_agmg+1;
    ctrl_loc=ctrl_solver;
    % passing index_agmg+1 we will use one copy of agmg, allowing preprocess 
    ctrl_loc.init(ctrl_inner22.approach,...
		  ctrl_inner22.tolerance,...
		  ctrl_inner22.itermax,...
		  ctrl_inner22.omega,...
		  ctrl_inner22.verbose,...
		  'C+B1diag(A)^{-1}B2',...
		  index_agmg);
    
    inv_SCA.init(approx_SCA+controls.relax4inv22*speye(Nr,Nr),ctrl_loc);
    inv_SCA.info_inverse.label='schur_ca';
    inv_SCA.cumulative_iter=0;
    inv_SCA.cumulative_cpu=0;
    
    inverse_cpu=inv_SCA.init_cpu;
    inverse_block22 = @(x) -inv_SCA.apply(x);

  elseif(strcmp(approach_schurCA,'iterative'))
				% define matrix-vector operator
    SCA = @(x) (C*x + B2*(inv_A(B1T*x)));
    % assembly inverse of  SAC iterative solver with no precodnitioner
    inverse_block22 = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);
    inner_nequ=Nr;

  elseif(strcmp(approach_schurCA,'full'))
    tildeS=form_minus_SchurCA(inv_A,B1T,B2,C);

    inv_SCA=sparse_inverse;

    % we need a new seed for agmg, in case is used
    index_agmg=index_agmg+1;
    ctrl_loc=ctrl_sovler;
    
    ctrl_loc=ctrl_solver;
    % passing index_agmg+1 we will use one copy of agmg, allowing preprocess 
    ctrl_loc.init(ctrl_inner22.approach,...
		  ctrl_inner22.tolerance,...
		  ctrl_inner22.itermax,...
		  ctrl_inner22.omega,...
		  ctrl_inner22.verbose,...
		  'schur_ca',...
		  index_agmg);
    
    inv_SCA.init(tildeS,ctrl_loc);
    inv_SCA.cumulative_iter=0;
    inv_SCA.cumulative_cpu=0;
    
    inverse_cpu=inv_SCA.init_cpu;
    inverse_block22 = @(x) -inv_SCA.apply(x);
  end

  if (0)%controls.compute_eigen)

				% create (~S)^{-1} S
    fprintf('building (~S)^{-1} S')
    fprintf('aaaaaa')
    precS=zeros(Nr,Nr);
    for i=1:Nr
      if (mod(i,100))
	fprintf('\b\b\b\b\b\b%4.2f ',i/Nr*100)
      end
      precS(:,i)=-inverse_block22(dense_S(:,i));
    end

    eigenvalues=eig(precS);
    real_part=real(eigenvalues);
    imag_part=imag(eigenvalues);
    abs_eig=abs(eigenvalues);
    fprintf('real part eigen (~S)^{-1} S : min=%1.4e max=%1.4e\n', min(real_part),max(real_part))
    fprintf('imag part eigen (~S)^{-1} S : min=%1.4e max=%1.4e\n', min(imag_part),max(imag_part))
    fprintf('abs       eigen (~S)^{-1} S : min=%1.4e max=%1.4e\n', min(abs_eig),max(abs_eig))

    figure
    plot(eigenvalues,'o')
  end

  

  
  
  % Define action of preconditoner
  prec = @(x) SchurCA_based_preconditioner(x, inv_A,...
					   inverse_block22,...
					   @(y) B1T*y,...
					   @(z) B2*z,...
					   Np,Nr,...
					   controls.outer_prec,ncellphi);

  % solve
  jacobian = [A B1T; B2 -C];
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) jacobian*x, rhs, ctrl_outer, prec );
  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);  
  end
  
  outer_cpu=toc(outer_timing);

  if(strcmp(approach_schurCA,'iterative'))
    if (strcmp(approach_inverse_A,'full'))
      total_Ainversion=inverseA.cumulative_application;
    elseif (strcmp(approach_inverse_A,'block'))
      total_Ainversion=diag_block_invA(1).cumulative_application;
    end
    fprintf('total number of A^{-1} application=%d\n',total_Ainversion)
  end
   
  % get info
  inner_iter=0;%inv_SCA.cumulative_iter;
  outer_iter=uint32(info_J.iter);
  
  
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  if (0)
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

  end
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  if (strcmp(approach_inverse_A,'full'))
    inverseA.kill();
  elseif (strcmp(approach_inverse_A,'block'))
    for i=1:Nt
      diag_block_invA(i).kill();
    end
  end

  inv_SCA.kill();
  
elseif sol==13
  verbose=controls.verbose;
  % approach from : "A PRECONDITIONER FOR GENERALIZED SADDLE POINT PROBLEMS"
  % solve M=(A B1T) = H + S
  %         (-B2 C) 
  % where
  % H=(A  ) S=(     B1T)
  %   (  C)   (-B2     )
  % Krylov method for full system
  % + 
  % use P = 1/(2a)*(H+aI)(S+aI)  
  % 
  
  % cpu time preprocessing main solver
  tic;

  % set rows with maximal entries on A diagonal 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  % get vectors p and q and alpha such that p^T x + q^y = alpha
  % we build N vectors
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls,controls.manipulation_approach);
  if (swap_sign)
    vectors_x=-vectors_x;
    vectors_y=-vectors_y;
    alphas=-alphas;
  end
  % manipulate A,B1T and rhs with p,q,alpha
  if (controls.manipulate>0)
    [A,B1T,rhs]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas);
  end

 
  				% grounding
  if (indc>0)
    disp('GROUNDING')
    inode=indeces_global(1);
    [A,B1T,rhs]=grounding(A,B1T,rhs,inode,0);
  elseif(indc<0)
    irow=indeces_global(1);
    A(irow,1:ncellphi) = -vectors_x(:,1)';
    B1T(irow,:)=sparse(1,Nr);
    rhs(irow,:)=0.0;
  end

  if (controls.diagonal_scaling)
    diagA_scaling = sparse(1:Np,1:Np,(1.0./sqrt(spdiags(A,0)))',Np,Np);
    diagC_scaling=sparse(1:Nr,1:Nr,(1.0./sqrt(spdiags(C,0)))',Nr,Nr);
    size(diagA_scaling);
    A=diagA_scaling*A*diagA_scaling;
    B1T=diagA_scaling*B1T*diagC_scaling;
    B2=diagC_scaling*B2*diagA_scaling;
    C=diagC_scaling*C*diagC_scaling;

    rhs(1:Np)      =diagA_scaling*rhs(1:Np);
    rhs(1+Np:Np+Nr)=diagC_scaling*rhs(1+Np:Np+Nr);   
  end

  
  % change sign of second row
  matrixM = [A B1T; -B2 C];
  rhs(Np+1:Np+Nr)=-rhs(Np+1:Np+Nr);

  %get relazation paramter
  alpha=controls.alpha;
 
  % define inverse of H+alphaI
  %inverse of A
  approach_inverse_A=controls.approach_inverse_A;
  if ( strcmp(approach_inverse_A,'full'))
    invAalpha=sparse_inverse;
    invAalpha.init(A+alpha*speye(Np,Np),controls.ctrl_innerA);
    invAalpha.cumulative_iter=0;
    invAalpha.cumulative_cpu=0;
			
    inv_Aalpha = @(x) invAalpha.apply(x);
  elseif( strcmp(approach_inverse_A,'block'))
    % define array of sparse inverse
    diag_block_invA(Nt,1)=sparse_inverse;

    % create blocks
    nAi=ncellphi;
    ctrl_loc=ctrl_solver;
    for i=1:Nt
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(controls.ctrl_innerA.approach,...
		    controls.ctrl_innerA.tolerance,...
		    controls.ctrl_innerA.itermax,...
		    controls.ctrl_innerA.omega,...
		    controls.ctrl_innerA.verbose,...
		    sprintf('%s A%d',controls.ctrl_innerA.label,i));
      diag_block_invA(i).name=sprintf('inverse A%d',i);
      
      % get block add relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)+alpha*speye(nAi,nAi);
      
      % define inverse
      diag_block_invA(i).init( matrixAi,ctrl_loc);     

    end
    % define function
    inv_Aalpha  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end

  % inverse of C+alpha*I
  inv_Calpha=sparse_inverse;
  inv_Calpha.init(C+alpha*speye(Nr,Nr),controls.ctrl_innerC);

  % concatenate
  inv_Halpha = @(x) [inv_Aalpha(x(1:Np)); inv_Calpha.apply(x(Np+1:Np+Nr))];

  % define inverse of S+alphaI
  % We can explicitely form the matrix
  if (strcmp(controls.approach_inverse_S,'SAC'))
    S11=alpha*speye(Np,Np)+1/alpha*B1T*B2;
    inv_S11=sparse_inverse;
    inv_S11.init(S11,controls.ctrl_innerS);
    
    inv_Salpha = @(x) SchurAC_based_preconditioner(x, @(xin) inv_S11.apply(xin),...
						   @(y) -y/alpha ,B1T,-B2,...
						   'full',ncellphi);
  elseif (strcmp(controls.approach_inverse_S,'SCA'))
    S11=alpha*speye(Nr,Nr)+1/alpha*B2*B1T;
    inv_S11=sparse_inverse;
    inv_S11.init(S11,controls.ctrl_innerS);
				% Define action of preconditoner
    inv_Salpha = @(x) SchurCA_based_preconditioner(x, @(xin) xin/alpha,...
						   @(yin) inv_S11.apply(yin),...
						   @(y) B1T*y,...
						   @(z) -B2*z,...
						   Np,Nr,...
						   'full',ncellphi);
  end


				% define preconditioner
  if (strcmp(controls.approach_prec,'HS'))
    prec = @(x) inv_Salpha(inv_Halpha(x))/(2*alpha);
  elseif ( strcmp(controls.approach_prec,'SH'))
    prec = @(x) inv_Halpha(inv_Salpha(x))/(2*alpha);
  end
 
  %preprocess finished
  preprocess_cpu=toc;


  % solve 
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) matrixM*x, rhs, controls.ctrl_outer, prec,[],controls.left_right );
  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);  
  end
  
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

%print_info_solver(sol_stat)
%return
