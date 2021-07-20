function [d, sol_stat, normd,resume_msg] = solvesys(JF,F,controls)

% solve the linear system d = -JF\F, reducing the system by eliminating the
% additional variable s and imposing the condition d(indc) = 0 to fix the 
% constant of the potential dp


sol=controls.sol;
indc=controls.indc;

resume_msg='';

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


inner_iter2=0;
inner_iter3=0;

inner_nequ2=0;
inner_nequ3=0;

flag=0;
iter=0;
res=1e30;


% rhs
f=-F.p;
g=-F.r;
h=-F.s;

% Solve with the Schur complement A|C, using Matlab's backslash (or agmg)
swap_sign=controls.swap_sign;

				% reduced system to phi rho variables
A = sparse(JF.pp); B1T = sparse(JF.pr); B2 = sparse(JF.rp);
R = sparse(JF.rr); M   = sparse(JF.rs);
Ds = sparse(JF.sr); Dr = sparse(JF.ss);
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
  
  R=-R;
  M=-M;
  Ds=-Ds;
  Dr=-Dr;

  
  f=-f;
  g=-g;
  h=-h;
  
  f1=-f1;
  f2=-f2;
  
end

% assembly rhs
rhs=[f1;f2];
rhs_full=[f;g;h];

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
  verbose=controls.verbose;
			% fullyreduced system
  prep=tic;

  time_manipulate=tic;
  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
  
   % grounding ,manipualtion etc
  [A,B1T,B2,rhs,B1T_perturbation]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);

  
% scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs

  if ( controls.permute)
    permutation = symrcm(A(1:ncellphi,1:ncellphi));
    p=permutation;
    for i=1:Nt
      range=[(i-1)*ncellphi+1:i*ncellphi];
      A(range,range) = A((i-1)*(ncellphi)+p,(i-1)*ncellphi+p);
      B1T(range,:)   = B1T((i-1)*(ncellphi)+p,:);
      B2(:,range)    = B2(:,(i-1)*(ncellphi)+p);
      rhs(range)     = rhs((i-1)*(ncellphi)+p);
    end
  end


  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end


  
  cpu_manipulate=toc(time_manipulate);
  msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
  if (verbose>1)
    fprintf('%s\n',msg);
    fprintf(controls.logID,'%s\n',msg);
  end


  
  time_invC_builS=tic;
				% init diag(C)^{-1}
  inverseC_approach=controls.inverseC_approach;
  if ( inverseC_approach==1)
    
    invC = sparse_inverse;

    ctrl_loc=ctrl_solver;
    index_agmg=1;
    ctrl_loc.init(controls.ctrl_inner22.approach,...
		  controls.ctrl_inner22.tolerance,...
		  controls.ctrl_inner22.itermax,...
		  controls.ctrl_inner22.omega,...
		  controls.ctrl_inner22.verbose,...
		  'C',index_agmg);
    
    invC.init(C,ctrl_loc);

    inverseC = @(x) invC.apply(x);

    diagC = spdiags(C,0);
    %if ( verbose >= 1)
    fprintf('(%9.4e <= C <= %18.12e) \n',min(diagC),max(diagC))
    %end
    invdiagC = sparse(1:Nr,1:Nr,(1.0./diagC)',Nr,Nr);

    
				% assembly schur AC S=A+B1
				% reduced to system only in phi
    
    S = A+B1T*(invdiagC*B2);  
    fp = rhs(1:Np)+B1T*(inverseC(rhs(Np+1:Np+Nr)));
    
  elseif ( inverseC_approach==2)
				%  C=-(R-Dr^{-1}M) Ds
				% ~C=-(D(rho)R - MD(S))
				% thus
				% C=Dr^{-1} ( ~C )
    tildeC=-Dr*R+M*Ds;

    inv_tildeC = sparse_inverse;

    ctrl_loc=ctrl_solver;
    index_agmg=1;
    ctrl_loc.init(ctrl_inner22.approach,...
		  ctrl_inner22.tolerance,...
		  ctrl_inner22.itermax,...
		  ctrl_inner22.omega,...
		  ctrl_inner22.verbose,...
		  'C',index_agmg);
    
    inv_tildeC.init(tildeC,ctrl_loc);

    
    
				% C^{-1}=(~C)^{-1})*Dr
    inverseC = @(x) inv_tildeC.apply(Dr*x);

    diagCtilde = spdiags(tildeC,0);
    if ( verbose >= 1)
      fprintf('(%9.4e <= Dr <= %18.12e) \n',min(spdiags(Dr,0)),max(spdiags(Dr,0)))
      fprintf('(%9.4e <= Ds <= %18.12e) \n',min(spdiags(Ds,0)),max(spdiags(Ds,0)))
      fprintf('(%9.4e <= diag(~C) <= %18.12e) \n',min(diagCtilde),max(diagCtilde))
    end
    inv_diag_tildeC = sparse(1:Nr,1:Nr,(1.0./diagCtilde)',Nr,Nr);

    S = A+B1T*((inv_diag_tildeC*Dr)*B2);     
  end

  cpu_invC_builS=toc(time_invC_builS);
  msg=sprintf('invC + BUILD S=%1.4e',cpu_invC_builS);
  if ( verbose >= 2)
    fprintf('%s\n',msg);
  end
  fprintf(controls.logID,'%s\n',msg);

  
  if ( verbose >= 2)
    fprintf('INIT S=A+B1T C^{-1} B2\n')
  end
  

 % assembly (approcimate) inverse of SAC=A+B1T diag(C)^{-1} B2 
  relax4prec=controls.relax4inv11;
  relax4prec=0.0;
  debug=0;

  

  build_S=tic;
  if ( verbose >= 2)
    fprintf('INIT inverse S\n')
  end
  inv_S=sparse_inverse;
  ctrl_loc=ctrl_solver;
  index_agmg=index_agmg+1;
  ctrl_loc.init(controls.ctrl_inner11.approach,...
		controls.ctrl_inner11.tolerance,...
		controls.ctrl_inner11.itermax,...
		controls.ctrl_inner11.omega,...
		controls.ctrl_inner11.verbose,...
		'SAC',index_agmg,...
		controls.logID);
  inv_S.init(S+relax4prec*speye(Np,Np),ctrl_loc);
  inv_S.cumulative_iter=0;
  inv_S.cumulative_cpu=0;
  inv_S.dimblock=ncellphi;
  inverse_cpu=inv_S.init_cpu;

  preprocess_S=toc(build_S);
  prec_cpu=preprocess_S;

  preprocess_cpu=toc(prep);
  
 % build rhs reduced
  fp = rhs(1:Np)+B1T*(inverseC(rhs(Np+1:Np+Nr)));


  % apply inverse S
  linsys=tic;
  dp=inv_S.apply(fp);

  
  % get y=C^{-1}(B2*x-g)
  dr=inverseC(B2*dp-rhs(Np+1:Np+Nr));
  d=[dp;dr];
  
  

  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);
  end


  if (controls.permute)
    for i=1:Nt
      range=[(i-1)*ncellphi+1:i*ncellphi];
      d((i-1)*(ncellphi)+p) = d(range);
    end
  end

  
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];

  outer_cpu = toc(linsys);


  
  flag=inv_S.info_inverse.flag;
  relres=inv_S.info_inverse.realres;
  iter=inv_S.info_inverse.iter;

  outer_iter=iter;
  
  inner_iter=0;
  inner_nequ=0;

  
  normd = norm(d);

  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.2e | assembly=%1.2e invS=%1.2e',...
		     flag,relres,resp,resr,ress, inv_S.cumulative_iter,outer_cpu,...
		     preprocess_cpu,preprocess_S);


  inv_S.kill();
    

elseif sol==10
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
  % use P=( SAC B1T ) as preconditoner 
  %       ( 0    -C )

  % cpu time preprocessing main solver
  tic;

  time_manipulate=tic;
  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
 
  % grounding ,manipualtion etc
  [A,B1T,B2,rhs,B1T_perturbation]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);

  size(  B1T_perturbation)
  
  % scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end

  
    
  cpu_manipulate=toc(time_manipulate);
  msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
  if (verbose>1)
    fprintf('%s\n',msg);
    fprintf(controls.logID,'%s\n',msg);
  end


  
  time_invC_builS=tic;
  % init diag(C)^{-1}
  inverseC_approach=controls.inverseC_approach;
  if ( inverseC_approach==1)
    
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

    inverseC = @(x) invC.apply(x);

    diagC = spdiags(C,0);
    if ( verbose >= 1)
      fprintf('(%9.4e <= C <= %18.12e) \n',min(diagC),max(diagC))
    end
    invdiagC = sparse(1:Nr,1:Nr,(1.0./diagC)',Nr,Nr);

   
    % assembly schur AC S=A+B1
    % reduced to system only in phi
    
    S = A+B1T*(invdiagC*B2);  
    fp = rhs(1:Np)+B1T*(inverseC(rhs(Np+1:Np+Nr)));
    
  elseif ( inverseC_approach==2)
    %  C=-(R-Dr^{-1}M) Ds
    % ~C=-(D(rho)R - MD(S))
    % thus
    % C=Dr^{-1} ( ~C )
    tildeC=-Dr*R+M*Ds;

    inv_tildeC = sparse_inverse;

    ctrl_loc=ctrl_solver;
    index_agmg=1;
    ctrl_loc.init(ctrl_inner22.approach,...
		  ctrl_inner22.tolerance,...
		  ctrl_inner22.itermax,...
		  ctrl_inner22.omega,...
		  ctrl_inner22.verbose,...
		  'C',index_agmg);
    
    inv_tildeC.init(tildeC,ctrl_loc);

				
    	
    % C^{-1}=(~C)^{-1})*Dr
    inverseC = @(x) inv_tildeC.apply(Dr*x);

    diagCtilde = spdiags(tildeC,0);
    if ( verbose >= 1)
      fprintf('(%9.4e <= Dr <= %18.12e) \n',min(spdiags(Dr,0)),max(spdiags(Dr,0)))
      fprintf('(%9.4e <= Ds <= %18.12e) \n',min(spdiags(Ds,0)),max(spdiags(Ds,0)))
      fprintf('(%9.4e <= diag(~C) <= %18.12e) \n',min(diagCtilde),max(diagCtilde))
    end
    inv_diag_tildeC = sparse(1:Nr,1:Nr,(1.0./diagCtilde)',Nr,Nr);

    S = A+B1T*((inv_diag_tildeC*Dr)*B2);     
  end

  cpu_invC_builS=toc(time_invC_builS);

  
  
  figure
  spy(S)


  ps = symrcm(S);

  S=S(ps,ps);

  figure
  spy(S)

  plot(ps)

  
  
  
  return
  
  msg=sprintf('invC + BUILD S=%1.4e',cpu_invC_builS);
  if ( verbose >= 2)
    fprintf('%s\n',msg);
  end
  fprintf(controls.logID,'%s\n',msg);

  
  if ( compute_eigen)
    eigvalues=eig(full(S));
    fprintf('eig(S) min %8.4e max %8.4e  \n', min(eigvalues), max(eigvalues))
  end
  
  % grounding the solution
  if ( verbose >= 2)
    fprintf('INIT S=A+B1T C^{-1} B2\n')
  end
 

  % assembly (approcimate) inverse of SAC=A+B1T diag(C)^{-1} B2 
  approach=controls.extra_info;  
  relax4prec=controls.relax4inv11;
  relax4prec=0.0;
  debug=0;

  timeS=tic;
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
		  'SAC',index_agmg,...
		  controls.logID);
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
    sol_rho=inverseC(B2*sol_phi-rhs(Np+1:Np+Nr));
    fprintf(' res |S x -f_r|/f_r = %1.4e \n',norm(S*sol_phi-fp)/norm(fp));
    
    
    jacobian=[A B1T;B2 -C];
    d=[sol_phi;sol_rho];
    
    
    fprintf(' res modified = %1.4e \n',norm(jacobian*d-rhs)/norm(rhs))
  end
  preprocess_S=toc(timeS);

  prec_cpu=cpu_invC_builS+preprocess_S;

  
  msg=sprintf('BUILD invS = %1.4e', preprocess_S);
  if ( verbose >= 2)
    fprintf('%s\n',msg);
  end
  fprintf(controls.logID,'%s\n',msg);

  preprocess_cpu=toc;
  
  %
  % Define action of preconditoner
  % 
  prec = @(x) SchurAC_based_preconditioner(x, inverse_block11, inverseC ,...
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

  if (0)    
    test_vectors(d,vectors_x,vectors_y,alphas,N,Np,Nr);
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

  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);

    %rhs(1:Np)      =diagA_scaling\rhs(1:Np);
    %rhs(1+Np:Np+Nr)=diagC_scaling\rhs(1+Np:Np+Nr);
  end
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];

  norm(jacobian(1:Np,:)*d(1:Np+Nr)-rhs(1:Np))/norm(rhs(1:Np))

  
  normd = norm(d);

  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.2e | nequ=%d - inner iter=%d | assembly=%1.2e invS=%1.2e',...
	  info_J.flag,relres,resp,resr,ress,outer_iter,outer_cpu,...
	  inner_nequ, inner_iter, ...
	  preprocess_cpu,preprocess_S);

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

  time_manipulate=tic;
  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
 
  % grounding ,manipualtion etc
  [A,B1T,B2,rhs]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);

  % scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end
  
  cpu_manipulate=toc(time_manipulate);
  msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
  if (verbose>1)
    fprintf('%s\n',msg);
    fprintf(controls.logID,'%s\n',msg);
  end


  time_prec=tic;
  % define approxiamte inverse of (~A)^{-1}  
  inverseA_approach=controls.extra_info;
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    invA=sparse_inverse;
    invA.init(A+controls.relax4inv11*speye(Np,Np),ctrl_inner11);
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
    fprintf('APPROACH %s \n',approach)
  end
  
  build_S=tic;
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

  cpu_assembly_inverseS=toc(build_S);

  prec_cpu=toc(time_prec);

  
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

  if (1)    
    test_vectors(d,vectors_x,vectors_y,alphas,N,Np,Nr);
  end

  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| S : nequ=%d inn=%d ass=%1.2e',...
   		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
   		     inv_SCA.nequ, inv_SCA.cumulative_iter,cpu_assembly_inverseS);

  inv_SCA.kill();

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

  time_manipulate=tic;
  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
 
  % grounding ,manipualtion etc
  [A,B1T,B2,rhs,B1T_perturbation]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);

  % scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end
  
  cpu_manipulate=toc(time_manipulate);
  msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
  if (verbose>1)
    fprintf('%s\n',msg);
    fprintf(controls.logID,'%s\n',msg);
  end

  time_prec=tic;


  % define inverse of (~A)^{-1}
  relax4prec=controls.relax4inv11;
  index_agmg=0;
  
  % inverse A action
  build_inverseA=tic;
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

    inverse11  = @(x) inverseA.apply(x);
    
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
		    controls.ctrl_inner11.verbose*(i==1),...
		    sprintf('%s A%d',ctrl_inner11.label,i),...
		    index_agmg);
      diag_block_invA(i).name=sprintf('inverse A%d',i);
      
      % get block add relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)+relax4prec*speye(nAi,nAi);
      
				% define inverse
      diag_block_invA(i).init( matrixAi,ctrl_loc);     

    end
% define function A shifter( inverse A grounded ( rhs set to zero in some nodes )
    inverse11  = @(x) apply_inverseAtilde(diag_block_invA,A,vectors_x,indeces_global,indeces_local,x);
  end

  cpu_assembly_inverseA=toc(build_inverseA);
  preprocess_cpu=toc;

  %
				%
  build_inverseS=tic;
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
    SCA = @(x) (C*x + B2*(inverse11(B1T*x)));
    % assembly inverse of  SAC iterative solver with no precodnitioner
    inverse_block22 = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);
    inner_nequ=Nr;

  elseif(strcmp(approach_schurCA,'iterative+SwithdiagA'))
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
    ctrl_loc.init(controls.ctrl_inner22inner.approach,...
		  controls.ctrl_inner22inner.tolerance,...
		  controls.ctrl_inner22inner.itermax,...
		  controls.ctrl_inner22inner.omega,...
		  controls.ctrl_inner22inner.verbose,...
		  'C+B1diag(A)^{-1}B2',...
		  index_agmg);
    
    inv_SCA.init(approx_SCA+controls.relax4inv22*speye(Nr,Nr),ctrl_loc);
    inv_SCA.info_inverse.label='schur_ca';
    inv_SCA.cumulative_iter=0;
    inv_SCA.cumulative_cpu=0;
    
    
				% define matrix-vector operator
    SCA = @(x) (C*x + B2*(inv_A(B1T*x)));
    % assembly inverse of  SAC iterative solver with no precodnitioner
    inverse_block22 = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) inv_SCA.apply(z));
    inner_nequ=Nr;
  elseif(strcmp(approach_schurCA,'timelaplacian'))
			       % define explicitely inverse of diag(A)

    deltat=1/(N+1);

	    % coarse A with PAPT ( maybe only diag(block A) is enough)

   
 

   
    
    


    DiagA = spdiags(A,0);
    invDiagA = sparse(1:Np,1:Np,(1.0./DiagA)',Np,Np);
	
				% define time laplacian



    
    B1T_time=-JF.B1T_time;


    B1T_space=-JF.B1T_space;
    DtT=-JF.DtT;
    if (~controls.indc==0)
      irow=indeces_global(1);
      B1T_time(irow,:)=sparse(1,Nr);
      B1T_space(irow,:)=sparse(1,Nr);
      DtT(irow,:)=sparse(1,Nr);
    end

    P=JF.It(1:Np,1:Nr);
    coarseA=(P')*JOC.Mxt*A*JOC.Mxt*(P)/3.0;

    
    time_laplacian=(DtT')*DtT;
    original=(B1T_time')*(A\B1T_time);
    size(coarseA)
    size(original)
    my=coarseA\time_laplacian;

    time_laplacian=(B1T_time')*(A\B1T_time);
    %figure
    %imagesc(original);

    %figure
    %imagesc(my);

    
    test=full(my-original);
    norm(test)
    eigenvalues=study_eigenvalues(test, '(commutator)',1);
    
    figure
    plot(eigenvalues,'o')
    return
    
    time_space=(B1T_space')*(A\(B1T_time ));
    %space=(JF.B1T_space')*(A\JF.B1T_space);
    space=(B1T_space')*(A\(B1T_space))+time_space+time_space';
    
    if (controls.diagonal_scaling)
      time_laplacian=diagC_scaling*time_laplacian*diagC_scaling;
      space=diagC_scaling*space*diagC_scaling;
    end

    SCA_left_right='left';
    if (strcmp(SCA_left_right,'right'))
      
      approx_SCA = C*coarseA+time_laplacian+space*coarseA;

			% we need a new seed for agmg, in case is used
      index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
% passing index_agmg+1 we will use one copy of agmg, allowing preprocess
      ctrl_loc.init(controls.ctrl_inner22.approach,...
		    controls.ctrl_inner22.tolerance,...
		    controls.ctrl_inner22.itermax,...
		    controls.ctrl_inner22.omega,...
		    controls.ctrl_inner22.verbose,...
		    'timelaplacian',...
		    index_agmg);

				% assembly approximate inverse
      inv_SCA=sparse_inverse;
      inv_SCA.init(approx_SCA+controls.relax4inv22*speye(Nr,Nr),ctrl_loc);
      inv_SCA.info_inverse.label='schur_ca';
      inv_SCA.cumulative_iter=0;
      inv_SCA.cumulative_cpu=0;

      inverse_block22 = @(x) -coarseA*inv_SCA.apply(x);

    elseif (strcmp(SCA_left_right,'left'))
      
      approx_SCA = coarseA*C+time_laplacian+coarseA*space;

			% we need a new seed for agmg, in case is used
      index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
% passing index_agmg+1 we will use one copy of agmg, allowing preprocess
      ctrl_loc.init(controls.ctrl_inner22.approach,...
		    controls.ctrl_inner22.tolerance,...
		    controls.ctrl_inner22.itermax,...
		    controls.ctrl_inner22.omega,...
		    controls.ctrl_inner22.verbose,...
		    'timelaplacian',...
		    index_agmg);

				% assembly approximate inverse
      inv_SCA=sparse_inverse;
      inv_SCA.init(approx_SCA+controls.relax4inv22*speye(Nr,Nr),ctrl_loc);
      inv_SCA.info_inverse.label='schur_ca';
      inv_SCA.cumulative_iter=0;
      inv_SCA.cumulative_cpu=0;

      inverse_block22 = @(x) -inv_SCA.apply(coarseA*x);
    end
     
      
    if (1)
      %space=sparse((JF.B1T_space')*(A\JF.B1T_space));
      dense_S=(C+B2*(A\B1T));

      

      %invDiagS = sparse(1:Nr,1:Nr,(1.0./sqrt(spdiags(dense_S,0)))',Nr,Nr);
      %spacediag=sqrt(spdiags(space,0));
      %Dspace = sparse(1:Nr,1:Nr,(spacediag)',Nr,Nr);
      %invDspace = sparse(1:Nr,1:Nr,(1.0./spacediag)',Nr,Nr);
      

      %figure;
      %colorbar;
				%imagesc(invDiagS*dense_S*invDiagS);
      %scaled=invDspace*space*invDspace;
      
      %droped=scaled;
      %droped(abs(droped)<1e-2)=0.0;
      %fprintf('nnz per row %f , %d - sparsity=%1.4e\n', nnz(droped)/Nr ,Nr/N, (nnz(droped)/(Nr^2))*100)
      %imagesc(droped);
      %res=Dspace*(droped-scaled)*Dspace;

      %test=full(dense_S\(approx_SCA));
		   % Create figure
      test=prec_times_matrix(inverse_block22,-dense_S);
      eigenvalues=study_eigenvalues(test, '(S)^{-1}droped',1);

      figure
      plot(eigenvalues,'o')
      


      saveas(gcf,strcat(controls.basename,controls.sys_name,'denseS.png'));

      
      %out=prec_times_matrix(inverse_block22,-dense_S);
      
      %eigenvalues=study_eigenvalues(out, '(Candidate)^{-1}S',1);
      
      %plot(eigenvalues,'o')
      return
    end

	
  elseif(strcmp(approach_schurCA,'full'))
    out=prec_times_matrix(@(x) inverse11(x),B1T);
    tildeS=C+B2*out;

    diagtildeS=sparse(1:Nr,1:Nr,(1.0/sqrt(diag(tildeS)))',Nr,Nr)

    tildeS=diagtildeS*tildeS*diagtildeS;

    tildeS=sparse(tildeS(abs(tildeS)>1-3));
    
    spy(tildeS)

    nnz(tildeS)/Nr

    return
    

    

    inv_SCA=sparse_inverse;

    % we need a new seed for agmg, in case is used
    index_agmg=index_agmg+1;
    
    ctrl_loc=ctrl_solver;
% passing index_agmg+1 we will use one copy of agmg, allowing preprocess
    disp(ctrl_inner22.approach);
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
    inverse_block22 = @(x) -diagtildeS*inv_SCA.apply(x);
  end
  cpu_assembly_inverseS=toc(build_inverseS);


  prec_cpu=toc(time_prec);

  
  if (controls.compute_eigen)
			% we need a new seed for agmg, in case is used
    index_agmg=index_agmg+1;
    ctrl_loc=ctrl_solver;
    ctrl_loc.init(controls.ctrl_inner11.approach,...
		  controls.ctrl_inner11.tolerance,...
		  controls.ctrl_inner11.itermax,...
		  controls.ctrl_inner11.omega,...
		  controls.ctrl_inner11.verbose,...
		  'A',index_agmg);

    
    invA=sparse_inverse;
    invA.init(A,ctrl_loc);
    B1T_time=JF.B1T_time;
    B1T_space=JF.B1T_space;
    LT=JF.LT;

    B2_time=JF.B1T_time';
    B2_space=JF.B1T_space';
    
    if (controls.swap_sign)
      B1T_time=-B1T_time;
      B1T_space=-B1T_space;
      B2_time=-B2_time;
      B2_space=-B2_space;
      LT=-LT;
    end

    nullA=sparse(Np,Np);
    nullrhs=sparse(Np,1);
    nullB2=sparse(Nr,Np);
    if (~controls.indc==0)
      irow=indeces_global(1);
      B1T_time(irow,:)=sparse(1,Nr);
      B1T_space(irow,:)=sparse(1,Nr);
    end
    [mtt,mss,dense_S]=study_eigen_S(A,@(x) invA.apply(x),B2_time,B2_space,B1T_time,B1T_space,B1T_perturbation,C);%,JF.diagrho,LT,...
%				    strcat(controls.basename,controls.sys_name));

    return
    inv_SCA=sparse_inverse;

    % we need a new seed for agmg, in case is used
    index_agmg=index_agmg+1;
    
    ctrl_loc=ctrl_solver;
% passing index_agmg+1 we will use one copy of agmg, allowing preprocess
    disp(ctrl_inner22.approach);
    ctrl_loc.init(ctrl_inner22.approach,...
		  ctrl_inner22.tolerance,...
		  ctrl_inner22.itermax,...
		  ctrl_inner22.omega,...
		  ctrl_inner22.verbose,...
		  'schur_ca',...
		  index_agmg);
    
    inv_SCA.init(sparse(dense_S),ctrl_loc);
    inv_SCA.cumulative_iter=0;
    inv_SCA.cumulative_cpu=0;

    inverse_block22 = @(x) -inv_SCA.apply(x);
   
  end
  
 

  
  % Define action of preconditoner
  prec = @(x) SchurCA_based_preconditioner(x, inverse11,...
					   inverse_block22,...
					   @(y) B1T*y,...
					   @(z) B2*z,...
					   Np,Nr,...
					   controls.outer_prec,ncellphi);

  % solve
  jacobian = [A B1T; B2 -C];
  outer_timing=tic;

  info_J=info_solver;

  [d,info_J]=apply_iterative_solver(@(x) jacobian*x, rhs, ctrl_outer, prec );
  %d=jacobian\rhs;
  %figure
  %spy(jacobian)
  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);  
  end
  
  outer_cpu=toc(outer_timing);

  
 
  % get info
  inner_iter=0;%inv_SCA.cumulative_iter;
  outer_iter=uint32(info_J.iter);
  
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  if (0)    
    test_vectors(d,vectors_x,vectors_y,alphas,N,Np,Nr);
  end

  

  
  %if(strcmp(approach_schurCA,'iterative'))
  if (strcmp(approach_inverse_A,'full'))
    total_Ainner=inverseA.cumulative_iter;
    total_Ainversion=inverseA.cumulative_application;
    inner_nequ_A=inverseA.nequ;
    inverseA.kill();
  elseif (strcmp(approach_inverse_A,'block'))
    total_Ainversion=0;
    total_Ainner=0;
    inner_nequ_A=diag_block_invA(1).nequ;
    
    for i=1:Nt
      total_Ainner=total_Ainner+diag_block_invA(i).cumulative_iter;
      total_Ainversion=total_Ainversion+diag_block_invA(i).cumulative_application;
      diag_block_invA(i).kill();
    end
  end
  %fprintf('total number of A^{-1} application=%d\n',total_Ainversion)
  %end
  
  if (strcmp(approach_schurCA,'diagA'))
    inner_nequ_S=inv_SCA.nequ;
    total_Sinner=inv_SCA.cumulative_iter;
    total_Sinversion=inv_SCA.cumulative_application;
    inv_SCA.kill();
  elseif (strcmp(approach_schurCA,'iterative'))
    inner_nequ_S=Nr;
    inner_nequ_S=0;
    total_Sinner=0;
  elseif (strcmp(approach_schurCA,'timelaplacian'))
    inner_nequ_S=Nr;
    total_Sinner=inv_SCA.cumulative_iter;
  end

  

  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| %s A: nequ=%d inn.=%d ass=%1.2e | S : nequ=%d inn=%d ass=%1.2e',...
	  info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
	  approach_inverse_A,inner_nequ_A, total_Ainner,cpu_assembly_inverseA, ...
	  inner_nequ_S, total_Sinner,cpu_assembly_inverseS);

  inner_nequ=inner_nequ_A;
  inner_iter=total_Ainner;
  
  inner_nequ2=inner_nequ_S;
  inner_iter2=total_Sinner;
  
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

  time_manipulate=tic;
  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
 
  % grounding ,manipualtion etc
  [A,B1T,B2,rhs]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);

  % scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end
  
  cpu_manipulate=toc(time_manipulate);
  msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
  if (verbose>1)
    fprintf('%s\n',msg);
    fprintf(controls.logID,'%s\n',msg);
  end

  
  
  % change sign of second row
  matrixM = [A B1T; -B2 C];
  rhs(Np+1:Np+Nr)=-rhs(Np+1:Np+Nr);

  %get relazation paramter
  alpha=controls.alpha;

  time_prec=tic;
  timeA=tic;
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

  cpu_assembly_inverseA=toc(timeA);

  % inverse of C+alpha*I
  inv_Calpha=sparse_inverse;
  inv_Calpha.init(C+alpha*speye(Nr,Nr),controls.ctrl_innerC);

  % concatenate
  inv_Halpha = @(x) [inv_Aalpha(x(1:Np)); inv_Calpha.apply(x(Np+1:Np+Nr))];

  % define inverse of S+alphaI
  % We can explicitely form the matrix
  timeS=tic;
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
  cpu_assembly_inverseS=toc(timeS);
  

				% define preconditioner
  if (strcmp(controls.approach_prec,'HS'))
    prec = @(x) inv_Salpha(inv_Halpha(x))/(2*alpha);
  elseif ( strcmp(controls.approach_prec,'SH'))
    prec = @(x) inv_Halpha(inv_Salpha(x))/(2*alpha);
  end
 
  %preprocess finished
  preprocess_cpu=toc;

  prec_cpu=toc(time_prec);

  % solve 
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) matrixM*x, rhs, controls.ctrl_outer, prec,[],controls.left_right );
  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);  
  end
  
  outer_cpu=toc(outer_timing);

  if ( strcmp(approach_inverse_A,'full'))
    inner_iter_A=invAalpha.cumulative_iter;
    inner_nequ_A=invAalpha.nequ;
  elseif ( strcmp(approach_inverse_A,'block'))
    inner_iter_A=0;
    inner_nequ_A=diag_block_invA(i).nequ;
    for i=1:Nt
      inner_iter_A=inner_iter_A+diag_block_invA(i).cumulative_iter;
    end
  end

  inner_iter_S=inv_S11.cumulative_iter;
  inner_nequ_S=inv_S11.nequ;
  
  inner_iter=0;  
  
  flag=info_J.flag;
  relres=info_J.res;
  outer_iter=info_J.iter;


  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| %s A: nequ=%d inn.=%d avg=%d ass=%1.2e | C : %s inn.=%d | S : nequ=%d inn=%d ass=%1.2e',...
		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
		     approach_inverse_A,inner_nequ_A,inner_iter_A,round(inner_iter_A/(Nt*outer_iter)),cpu_assembly_inverseA, ...
		     inv_Calpha.ctrl.label,inv_Calpha.cumulative_iter,...
		     inner_nequ_S, inner_iter_S,cpu_assembly_inverseS);

  inner_nequ=inner_nequ_A;
  inner_iter=inner_iter_A;
  
  inner_nequ2=inner_nequ_S;
  inner_iter2=inner_iter_S;

  inner_nequ3=Nr;
  inner_iter3=inv_Calpha.cumulative_iter;

  
  
  % free memory
  if ( strcmp(approach_inverse_A,'full'))
    invAalpha.kill();
  elseif ( strcmp(approach_inverse_A,'block'))
    for i=1:Nt
      diag_block_invA(i).kill();
    end
  end
  inv_Calpha.kill();
  inv_S11.kill();
  
  
elseif sol==14
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
  preprocess=tic;

  time_manipulate=tic;
  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
 
  % grounding ,manipualtion etc
  [A,B1T,B2,rhs,B1T_perturbation]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);

  % scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end
  
  cpu_manipulate=toc(time_manipulate);
  msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
  if (verbose>1)
    fprintf('%s\n',msg);
    fprintf(controls.logID,'%s\n',msg);
  end

  gamma=controls.gamma;
				% init diag(C)^{-1}
  
  W_approach=controls.W_approach;
  if ( strcmp(W_approach,'C') )
    diagW=spdiags(C,0);
    W=sparse(1:Nr,1:Nr,diagW,Nr,Nr);
    invW=sparse(1:Nr,1:Nr,(1.0./diagW)',Nr,Nr);
  elseif ( strcmp(W_approach,'cutC') )
    diagW=spdiags(C,0);
    fprintf('(%9.4e <= w=diag(~C) <= %18.12e) \n',min(diagW),max(diagW))
    diagW(abs(diagW)<controls.lower_bound)=controls.lower_bound;
    diagW(abs(diagW)>controls.upper_bound)=controls.upper_bound;
    W=sparse(1:Nr,1:Nr,diagW,Nr,Nr);
    invW=sparse(1:Nr,1:Nr,(1.0./diagW)',Nr,Nr);
  elseif ( strcmp(W_approach,'Mass') )
    W=-JF.rs;
    invW=sparse(1:Nr,1:Nr,(1.0./spdiags(W,0))',Nr,Nr);
  elseif ( strcmp(W_approach,'MassgammaC') )
    W=-JF.rs+gamma*sparse(1:Nr,1:Nr,spdiags(C,0),Nr,Nr);
    invW=sparse(1:Nr,1:Nr,(1.0./spdiags(W,0))',Nr,Nr);
  end
    

  build_inv_augA=tic;
				% augemented system

  
  
  augA=A+gamma*B1T*invW*B2;
  augB1T=B1T*(speye(Nr,Nr)-gamma*invW*C);
  %norm(nonzeros(augB1T))
  aug_f=rhs(1:Np)+gamma*B1T*invW*rhs(Np+1:Np+Nr);
  augrhs=zeros(Np+Nr,1);
  augrhs(1:Np)=aug_f; %f_gamma
  augrhs(1+Np:Np+Nr)=rhs(Np+1:Np+Nr); %g

  % init inverse of augA

  tA=tic;
  ctrl_loc=ctrl_solver;
  index_agmg=0;
  index_agmg=index_agmg+1;
  ctrl_loc.init(ctrl_inner11.approach,...
		ctrl_inner11.tolerance,...
		ctrl_inner11.itermax,...
		ctrl_inner11.omega,...
		ctrl_inner11.verbose,...
		'augA',index_agmg);

  inv_augA=sparse_inverse;
  inv_augA.init(augA,ctrl_loc);
  inv_augA.cumulative_iter=0;
  inv_augA.cumulative_cpu=0;
  inv_augA.dimblock=ncellphi;
				%inv_S.info();
  
  inverse11 = @(x) inv_augA.apply(x);

  inner_nequ=Np;
  
  cpu_augA=toc(tA);

  
  tS=tic;
  % init diag(W)^{-1}
  S_approach=controls.S_approach;
  if ( strcmp(S_approach,'Mass') )
    ctrl_loc=ctrl_solver;
    index_agmg=index_agmg+1;
    ctrl_loc.init(ctrl_inner22.approach,...
		  ctrl_inner22.tolerance,...
		  ctrl_inner22.itermax,...
		  ctrl_inner22.omega,...
		  ctrl_inner22.verbose,...
		  'W',index_agmg);
    
    inv_S = sparse_inverse;
    inv_S.init(W,ctrl_loc);
   
    inverse22 = @(x) -gamma*inv_S.apply(x);

  elseif ( strcmp(S_approach,'GammaMassC') )
      ctrl_loc=ctrl_solver;
      index_agmg=index_agmg+1;
      ctrl_loc.init(ctrl_inner22.approach,...
		    ctrl_inner22.tolerance,...
		    ctrl_inner22.itermax,...
		    ctrl_inner22.omega,...
		    ctrl_inner22.verbose,...
		    'gammaMassC',index_agmg);
      
      inv_S = sparse_inverse;
      inv_S.init(1/gamma*W+C,ctrl_loc);
      inverse22 = @(x) -inv_S.apply(x);
  elseif ( strcmp(S_approach,'W') )
      ctrl_loc=ctrl_solver;
      index_agmg=index_agmg+1;
      ctrl_loc.init(ctrl_inner22.approach,...
		    ctrl_inner22.tolerance,...
		    ctrl_inner22.itermax,...
		    ctrl_inner22.omega,...
		    ctrl_inner22.verbose,...
		    'gammaMassC',index_agmg);
      
      inv_S = sparse_inverse;
      inv_S.init(1/gamma*W,ctrl_loc);
      inverse22 = @(x) -inv_S.apply(x);
  elseif ( strcmp(S_approach,'C') )
      ctrl_loc=ctrl_solver;
      index_agmg=index_agmg+1;
      ctrl_loc.init(ctrl_inner22.approach,...
		    ctrl_inner22.tolerance,...
		    ctrl_inner22.itermax,...
		    ctrl_inner22.omega,...
		    ctrl_inner22.verbose,...
		    'invC',index_agmg);
      
      inv_S = sparse_inverse;
      inv_S.init(C,ctrl_loc);
      inverse22 = @(x) -inv_S.apply(x);

  end
  cpu_S=toc(tS);

  prec_cpu=cpu_augA+ cpu_S;

  
  
  
  %
  % Define action of preconditoner
  %
  % assembly 2x2 prec
  prec = @(x) SchurCA_based_preconditioner(x, inverse11,...
					   inverse22,...
					   @(y) augB1T*y,...
					   @(z) B2*z,...
					   Np,Nr,...
					   controls.outer_prec,ncellphi);



  preprocess_cpu=toc(preprocess);

  
  % solve with fgmres 
  outer_timing=tic;
  augjacobian = [augA augB1T; B2 -C];
  if ( verbose >= 2)
    fprintf('START SOLVER \n')
  end
  [d,info_J]=apply_iterative_solver(@(x) augjacobian*x, augrhs, ctrl_outer, prec );
  if ( verbose >= 2)
    fprintf('END SOLVER \n')
  end
  info_J.print();
  outer_iter=uint64(info_J.iter);
  outer_cpu=toc(outer_timing);

  inner_iter=inv_augA.cumulative_iter;

  
  
  
  % store info
  flag=info_J.flag;
  iter=info_J.iter;

  
  
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);

  relres=ressys;

  S_approach=controls.S_approach;
  if ( strcmp(S_approach,'Mass') )
    inner_nequ2=inv_S.nequ;
    inner_iter2=inv_S.cumulative_iter;
  elseif ( strcmp(S_approach,'GammaMassC') )
    inner_nequ2=inv_S.nequ;
    inner_iter2=inv_S.cumulative_iter;
  elseif ( strcmp(S_approach,'C') )
    inner_nequ2=inv_S.nequ;
    inner_iter2=inv_S.cumulative_iter;

  end

  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.2e | Aug nequ=%d - iter=%d  built=%1.2e| S: nequ=%d - iter=%d  built=%1.2e',...
		     info_J.flag,relres,resp,resr,ress,outer_iter,outer_cpu,...
		     inner_nequ, inner_iter,cpu_augA, ...
		     inner_nequ2, inner_iter2,cpu_S);

  
  inner_nequ3=0;
  inner_iter3=0;


  % free memory (mandatory for agmg)
  inv_augA.kill();
  inv_S.kill();

elseif (sol==15)
			  % solve full system using recursive function

  % copy 
  copy_controls=controls;  
  copy_JF=JF;
  copy_F=F;

  copy_controls.sol=controls.block22_sol;
  
  % define globals preconditoner rhs and matrix vectors product
  pre_full = @(x) prec_solvesys(copy_JF,copy_F,copy_controls,x)

  if (controls.swap_sign)
    sign=-1
  else
    sign=1
  end
  rhs=sign*[-F.p;F.r;F.s];

  
  [d,info_J]=apply_iterative_solver(@(x) sign*J_times_vectors(copy_JF,x),rhs , copy_controls.ctrl_outer, prec_full )

elseif (sol==16)
  verbose=controls.verbose;
 % cpu time preprocessing main solver
  tic;

  time_manipulate=tic;
  % set rows 
  [indeces_global, indeces_local]=set_grounding_node(A,ncellphi); 
  [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
  
				% grounding ,manipualtion etc
  [A,B1T,B2,rhs,B1T_perturbation]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);

% scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end
  
  cpu_manipulate=toc(time_manipulate);
  msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
  if (verbose>1)
    fprintf('%s\n',msg);
    fprintf(controls.logID,'%s\n',msg);
  end

  % prec assembly
  time_assembly=tic;

  % inverse A
  build_inverseA=tic;

  index_agmg=0;
  % inverse A action
  if (strcmp(controls.approach_inverse_A,'full'))
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
    inverseA.init(A,ctrl_loc);
    inverseA.info_inverse.label='A';
    inverseA.cumulative_iter=0;
    inverseA.cumulative_cpu=0;

    inverse11  = @(x) inverseA.apply(x);
    
  elseif(strcmp(controls.approach_inverse_A,'block'))
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
		    sprintf('%s A%d',controls.ctrl_inner11.label,i),...
		    index_agmg);
      diag_block_invA(i).name=sprintf('inverse A%d',i);
      
      % get block add relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      
      % define inverse
      diag_block_invA(i).init( matrixAi,ctrl_loc);     

    end
% define function A shifter( inverse A grounded ( rhs set to zero in some nodes )
    inverse11  = @(x) apply_inverseAtilde(diag_block_invA,A,vectors_x,indeces_global,indeces_local,x);
  end

  cpu_assembly_inverseA=toc(build_inverseA);

  time_S=tic;
  % inverse S
  Q=JF.Mxt;
  Q=speye(Np,Np);
  invQ=sparse(1:Np,1:Np,(1./spdiags(Q,0))',Np,Np);
  
  % least square prec
  matrixF=(B2*invQ*B1T+controls.gamma*C);

  % init F^{-1}
  ctrl_loc=ctrl_solver;
  index_agmg=index_agmg+1;
  ctrl_loc.init(controls.ctrl_inner22.approach,...
		controls.ctrl_inner22.tolerance,...
		controls.ctrl_inner22.itermax,...
		controls.ctrl_inner22.omega,...
		controls.ctrl_inner22.verbose,...
		'F',index_agmg);

  inv_matrixF = sparse_inverse;
  inv_matrixF.init(matrixF,ctrl_loc);

  LSQ = @(x) inv_matrixF.apply(B2*invQ*A*invQ*B1T*(inv_matrixF.apply(x))) ;

  invdiagA=sparse(1:Np,1:Np,(1./spdiags(A,0))',Np,Np);
  Dmatrix=spdiags(B2*invdiagA*B1T+(1/controls.gamma)*C,0);

  inverse22 = @(x) -(LSQ(x) + controls.alpha*Dmatrix.*(x));

  if (controls.compute_eigen)
			% we need a new seed for agmg, in case is used
    index_agmg=index_agmg+1;
    ctrl_loc=ctrl_solver;
    ctrl_loc.init('direct',...
		  1e-14,...
		  1000,...
		  0,...
		  0,...
		  'A',index_agmg);

    accurate_inverseA = sparse_inverse;
    accurate_inverseA.init(A,ctrl_loc);
    
    out=prec_times_matrix(@(x) accurate_inverseA.apply(x), B1T);
    denseS=C+B2*out;

    %study_eigenvalues(denseS,'S',1);

    prec_matrix=prec_times_matrix(inverse22,-denseS);

    eigenvalues=study_eigenvalues(prec_matrix,'~S^{-1} S',1);%controls.logID);

    figure
    plot(eigenvalues,'o')
				
  
    saveas(gcf,strcat(controls.basename,'precS.png'));
    
    accurate_inverseA.kill();
  end

  

  cpu_assembly_inverseS=toc(time_S);
  
  % assembly 2x2 prec
  prec = @(x) SchurCA_based_preconditioner(x, inverse11,...
					   inverse22,...
					   @(y) B1T*y,...
					   @(z) B2*z,...
					   Np,Nr,...
					   controls.outer_prec,ncellphi);

  jacobian = [A B1T; B2 -C];
  cpu_assembly=toc(time_assembly);


  % solve
  info_J=info_solver;
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) jacobian*x, rhs, controls.ctrl_outer, prec,[],controls.left_right );
  if (controls.diagonal_scaling)
    d(1:Np)      =diagA_scaling*d(1:Np);
    d(1+Np:Np+Nr)=diagC_scaling*d(1+Np:Np+Nr);  
  end

  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  
  normd = norm(d);

  outer_iter=info_J.iter;
  outer_cpu=toc(outer_timing);



  
  % get info
  if ( strcmp(controls.approach_inverse_A,'full'))
    inner_iter_A=inverseA.cumulative_iter;
    inner_nequ_A=inverseA.nequ;
  elseif ( strcmp(controls.approach_inverse_A,'block'))
    inner_iter_A=0;
    inner_nequ_A=diag_block_invA(i).nequ;
    for i=1:Nt
      inner_iter_A=inner_iter_A+diag_block_invA(i).cumulative_iter;
    end
  end

  inner_iter_S=inv_matrixF.cumulative_iter;
  inner_nequ_S=inv_matrixF.nequ;

  

  
  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| %s A: nequ=%d inn.=%d ass=%1.2e | S : nequ=%d inn=%d ass=%1.2e',...
		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
		     controls.approach_inverse_A,inner_nequ_A, inner_iter_A,cpu_assembly_inverseA, ...
		     inner_nequ_S, inner_iter_S,cpu_assembly_inverseS);

  inner_iter=inner_iter_A;

  % free memory
  inv_matrixF.kill();
  if ( strcmp(controls.approach_inverse_A,'full'))
    inverseA.kill();
  elseif ( strcmp(controls.approach_inverse_A,'block'))
    for i=1:Nt
      diag_block_invA.kill();
    end
  end
  
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
		  'outer_cpu',  outer_cpu,...
		  'prec_cpu',prec_cpu,...
		  'inner_nequ2', inner_nequ2,...
		  'inner_iter2', inner_iter2,...
		  'inner_nequ3', inner_nequ3,...
		  'inner_iter3', inner_iter3);

%print_info_solver(sol_stat)
%return
