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

if (controls.swap_sign)
  sign=-1;
else
  sign=1;
end

				% reduced system to phi rho variables
A = sparse(JF.pp);
B1T = sparse(JF.pr-JF.ps*(JF.ss\JF.sr));
B2 = sparse(JF.rp);
R = sparse(JF.rr);
M   = sparse(JF.rs);
Ds = sparse(JF.sr);
Dr = sparse(JF.ss);
C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
temp=(JF.ss\h);
f1 = f-JF.ps*temp;
f2 = g-JF.rs*temp;

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


area2h=JF.area2h;

deltat=1/(N+1);
  
% build matrix W
W_mat=zeros(N,Nr);
for i = 1:N
  W_mat(i,1+(i-1)*ncellrho:i*ncellrho)=deltat*area2h';
end

%P=speye(Nr,Nr)-W_mat'*W_mat;

% define implicetly projector
factor=1/(deltat*norm(area2h))^2 
proj=@(rho) rho - factor * W_mat'*(W_mat*rho);

% e=zeros(Np,1);
% e(1:ncellphi)=1;
% temp=(B1T)'*e;
% % I/delta=I*Nt
% area2h=temp(1:ncellrho)*(1/Nt);
% figure
% plot(area2h)
% hold on
% plot(area2h1-area2h)
% hold off
% return

if sol==1
    
    % Solve using Matlab's backslash
    
    %eps = 1e-8;
    %eId = eps*speye(size(JF.pp));
  if (1)
    Sys = [A B1T; B2 -C];
    fp=f1;
    Sys(indc,:) = sparse(1,Np+Nr); Sys(:,indc) = sparse(Np+Nr,1); Sys(indc,indc) = 1; fp(indc)=0; 
    d = Sys\[fp;f2];

    res=norm(Sys*d-[f1;f2])/norm([f1;f2]);
    dp = d(1:Np); dr = d(Np+1:end);
    ds = Dr\(h-Ds*dr);
    d = [dp; dr; ds];

    normd = norm(d);
    resume_msg=sprintf('direct solver res=%1.2e',norm(res));
    
    
  else
    Sys = sparse(...
		  [JF.pp,        JF.pr, JF.ps;
		   JF.rp,        JF.rr, JF.rs;
		   sparse(Nr,Np), JF.sr, JF.ss]...
	    );
    
    rhs=-[f;g;h];
    Sys(indc,:) = sparse(1,Np+Nr+Nr); Sys(:,indc) = sparse(Np+Nr+Nr,1); Sys(indc,indc) = 1; rhs(indc)=0.0;
    d = Sys\rhs;
    res=norm(Sys*d-rhs)/norm(rhs);

    
    %eigenvalues= study_eigenvalues(full(Sys),'J_Full');
    %figure
    %plot(eigenvalues,'o')
    %return

				%dp = d(1:Np); dr = d(Np+1:end);
    %ds = Dr\(h-Ds*dr);
    %d = [dp; dr; ds];

  end

  if (0)
    % we enforse the rho increment to have zero mass

    
    % P=  (I_{m-1}
    %    -(a1 -a2 -a_{m-1})/a_m 
    mat_P=speye(ncellrho,ncellrho-1);   
    mat_P(ncellrho,1:ncellrho-1)=-area2h(1:ncellrho-1)'/area2h(ncellrho);

    % H=block diagonal (P) 
    mat_H = repmat({mat_P},1,N);
    mat_H = blkdiag(mat_H{:});


    s_B1T = B1T * mat_H;
    s_B2 = mat_H' * B2;
    s_C  = mat_H' * C * mat_H;

    s_rhs = [f1; mat_H' * f2];

    
    J_short=[A, s_B1T; s_B2, -s_C];

    d_short = J_short\s_rhs;

    dp = d_short(1:Np);
    dr = mat_H*d_short(Np+1:size(d_short,1));    
    ds = Dr\(h-Ds*dr);



    res_r=-(JF.rp * dp + JF.rr * dr + JF.rs * ds + F.r);
    dl=zeros(N,1);

    deltat=1/(N+1);
    factor=1/(deltat*norm(area2h))^2;
    for i = 1:N    
      dl(i)= factor * (area2h'* (res_r(1+(i-1)*ncellrho:i*ncellrho))) ;
      fprintf('dlambda(%d)=%1.2e\n',i,dl(i)) 
    end
    
    sum_dl=0;
    for i = 1:N
      sum_dl=sum_dl+dl(i);
      dp(1+(i-1)*ncellphi:i*ncellphi)=dp(1+(i-1)*ncellphi:i*ncellphi)+deltat*sum_dl;
    end
    
    d = [dp; dr; ds];

    [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
      
    
				% redefine rho dimension 

  end


  [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
  normd = norm(d);
  resume_msg = sprintf('direct solver res=%1.2e',resnorm);

  %figure
  %spy(Sys)
  %return
    
    inner_nequ=0;
    inner_iter=0;
    outer_iter=0;
    outer_cpu=0;
    prec_cpu=0;
    inner_nequ2;
    inner_iter2=0;
    inner_nequ3=0;
    inner_iter3=0;
    
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

  if ( controls.permute == 1)
    % define permutation
    perm_loc = symrcm(A(1:ncellphi,1:ncellphi));
    perm=(1:Np)';
    
    for i=1:Nt
      perm((i-1)*ncellphi+1:i*ncellphi)=(i-1)*(ncellphi)+perm_loc;
    end
    iperm(perm)=(1:Np)';
    
    % apply   fprintf('res=%1.2e\n',norm(S_old*x-b))

    A=A(perm,perm);
    B1T=B1T(perm,:);
    B2=B2(:,perm);
    rhs(1:N)=rhs(perm);
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
    %fprintf('(%9.4e <= C <= %18.12e) \n',min(diagC),max(diagC))
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
  S_old= S;

  
  figure
  spy(S)
  
  if ( controls.permute == 2)
    perm = symrcm(S);
    iperm(perm)=(1:Np)';
    % permute row and column
    S   = S(perm,perm);
    % permute row
    B1T = B1T(perm,:);
    % permute column
    B2 = B2(:,perm);
    % permute row
    rhs(1:Np) = rhs(perm);
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

  t=S_old*dp(iperm)-fp(iperm);
  fprintf('res=%1.2e\n',norm(t(iperm)))

  figure
  spy(S(iperm,iperm))

  

  b=fp(iperm);
  A = sign*sparse(JF.pp);
  B1T = sign*sparse(JF.pr-JF.ps*(JF.ss\JF.sr));
  B2 = sign*sparse(JF.rp);
  %S_old= A + B1T*(invdiagC*B2);
  norm(full(S_old-S(iperm,iperm)))
  x=dp(iperm);
  fprintf('res=%1.2e\n',norm(S_old*x-b))
  return
  

  if (controls.diagonal_scaling)
    dp = diagA_scaling*dp;
  end

  if (controls.diagonal_scaling)
    dr=diagC_scaling*dr;
  end


  % get y=C^{-1}(B2*x-g)
  dr = inverseC(B2*dp-rhs(Np+1:Np+Nr));

  
  if (controls.permute==1)
    dp = dp(perm);
  elseif (controls.permute==2)
    dp = dp(perm);
  end


  
  
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
    if ( strcmp(controls.assembly_S,'full') )
				%S = A+B1T*(invdiagC*B2);
      S = A+...
	  JF.B1T_time* (invdiagC*JF.B1T_time' )+...
	  JF.B1T_space*(invdiagC*JF.B1T_space')+...	 
	  JF.B1T_time* (invdiagC*JF.B1T_space')+...
	  JF.B1T_space*(invdiagC*JF.B1T_time');

    elseif ( strcmp(controls.assembly_S,'A_Mtt') ) 
      S = A+JF.B1T_time*(invdiagC*JF.B1T_time');
    elseif ( strcmp(controls.assembly_S,'A_Mtt_Mxx') ) 
      S = A+...
	  JF.B1T_time* (invdiagC*JF.B1T_time' )+...
	  JF.B1T_space*(invdiagC*JF.B1T_space');
    elseif ( strcmp(controls.assembly_S,'A_Mtx_Mxt') )
      S = A+...
	  JF.B1T_time* (invdiagC*JF.B1T_space')+...
	  JF.B1T_space*(invdiagC*JF.B1T_time');
    elseif ( strcmp(controls.assembly_S,'A_Mtt_Mtx_Mxt') )
      S = A+...
	  JF.B1T_time* (invdiagC*JF.B1T_time' )+...	  
	  JF.B1T_time* (invdiagC*JF.B1T_space')+...
	  JF.B1T_space*(invdiagC*JF.B1T_time');
    elseif ( strcmp(controls.assembly_S,'A_Mtt_Mtx_Mxt_blockdiagMxx') )
      block=sparse(Np,Np);
      temp=JF.B1T_space*(invdiagC*JF.B1T_space');
      for i=1:Nt
	block(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi)=...
	temp(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi);
      end
      S = A+...
	  JF.B1T_time* (invdiagC*JF.B1T_time' )+...	  
	  JF.B1T_time* (invdiagC*JF.B1T_space')+...
	  JF.B1T_space*(invdiagC*JF.B1T_time')+...
	  block;
    end 

    
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


  if ( strcmp(controls.cutS,'lower'))
    temp=JF.B1T_space*(invdiagC*JF.B1T_space');
    for i=1:N
      S(1+(i)*ncellphi:(i+1)*ncellphi,1+(i-1)*ncellphi:i*ncellphi)=...
      sparse(ncellphi,ncellphi);
    end
  end
    
  
  
  %figure
  %spy(S)


  %ps = symrcm(S);

  %S=S(ps,ps);

  %figure
  %spy(S)

  %plot(ps)

  
  
  
  %return
  
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
  if (0)
    time_manipulate=tic;
				% set rows 
    [indeces_global, indeces_local]=set_grounding_node(A,ncellphi);
    [vectors_x,vectors_y,alphas]=get_vectors_alphas_from_OC(JF,F,controls);
    
				% grounding ,manipualtion etc
    [A,B1T,B2,rhs]=preprocess_system(A,B1T,B2,rhs,indeces_global,vectors_x,vectors_y,alphas,controls);
    cpu_manipulate=toc(time_manipulate);
    msg=sprintf('MANIPULATE LIN. SYS= %1.4e',cpu_manipulate);
    if (verbose>1)
      fprintf('%s\n',msg);
      fprintf(controls.logID,'%s\n',msg);
    end
  end

  %[A,B1T,B2,rhs]=grounding(A,B1T,B2,rhs,controls.indc,0);
  
% scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end
  

  

  time_prec=tic;
  % define approxiamte inverse of (~A)^{-1}  
  inverseA_approach=controls.inverse11;
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    invA=sparse_inverse;
    invA.init(A+controls.relax4inv11*speye(Np,Np),ctrl_inner11);
    invA.cumulative_iter=0;
    invA.cumulative_cpu=0;

    % define function
    inv_A = @(x) invA.apply(x);
  elseif( strcmp(inverseA_approach,'block'))
    % partion matrix 
    nAi=ncellphi;
   
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
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + controls.relax4inv11*speye(nAi,nAi) ;
      diag_block_invA(i).init(matrixAi,ctrl_loc);     
    end

    above_diag_block_A= cell(N, 1);
    nAi=ncellphi;
    for i=1:N
      above_diag_block_A{i}=A((i-1)*nAi+1 : (i)*nAi , (i)*nAi+1 : (i+1)*nAi);
    end
    i=Nt;

				% define function
    inv_A = @(x) apply_block_triangular_inverse(diag_block_invA,above_diag_block_A,'U',x);
    %inv_A  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  elseif( strcmp(inverseA_approach,'diag'))
    % partion matrix 
    nAi=ncellphi;
   
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
    		    sprintf('BD%sA%d',ctrl_inner11.label,i));
      diag_block_invA(i).name=sprintf('inverse A%d',i);
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + controls.relax4inv11*speye(nAi,nAi) ;
      diag_block_invA(i).init(matrixAi,ctrl_loc);     
    end

    % define function
    inv_A  = @(x) projector(apply_block_diag_inverse(diag_block_invA,x),kernel);
  end

  
  
  % define explicitely inverse of diag(A)
  invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);
  

  if (controls.remove_imbalance)
    % we enforse the rho increment to have zero mass

    
    % P=  (I_{m-1}
    %    -(a1 -a2 -a_{m-1})/a_m 
    mat_P=speye(ncellrho,ncellrho-1);   
    mat_P(ncellrho,1:ncellrho-1)=-area2h(1:ncellrho-1)'/area2h(ncellrho);

    % H=block diagonal (P) 
    mat_H = repmat({mat_P},1,N);
    mat_H = blkdiag(mat_H{:});

    norm(mat_P'*area2h)
    figure

    B1T = B1T * mat_H;
    B2 = mat_H' * B2;
    C  = mat_H' * C * mat_H;

    % for i=1:Nt
    %   e=zeros(Np,1);
    %   e(1+(i-1)*ncellphi:i*ncellphi)=1;
    %   temp=B2*e;
    %   figure
    %   plot(temp)
    %   norm(temp)
    % end
    
    rhs = [rhs(1:Np); mat_H' * rhs(Np+1:Np+Nr)];

    % redefine rho dimension 
    Nr = (ncellrho-1)*N;
  end
  
  % assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
  % reduced to system only in rho
  SCA = (C+B2*invDiagA*B1T);

  % preprocess of precondtioner defienition
  preprocess_cpu=toc;
  
  debug=0;  
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
  if ( verbose >= 2)
    fprintf('DONE inverse S\n') 
  end 
  inverse_block22 = @(x) -inv_SCA.apply(x);

  cpu_assembly_inverseS=toc(build_S);

  %SCA = @(x) (C*x + B2*(invA(B1T*x)));
  %assembly inverse of  SAC iterative solver with no precodnitioner
  %inverse_block22 = @(y) -apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);
  

  % store time required to build prec
  prec_cpu=toc(time_prec);

  % set dimension of inner solver
  inner_nequ=Nr;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % END DUAL SCHUR COMPLEMENT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  
  % Define action of block preconditoner 
  prec = @(x) SchurCA_based_preconditioner(x, inv_A,...
					   inverse_block22,...
					   @(y) B1T*proj(y),...
					   @(z) B2*z,...
					   Np,Nr,...
					   controls.outer_prec,ncellphi);
  
  %proj = @(z) [projector(z(1:Np),kernel);z(Np+1:Np+Nr)]; 
  
 % solve
  if (0)
    % assembly full system
    jacobian = [A B1T; B2 -C];

    
    out=prec_times_matrix(prec,jacobian);
    norm(jacobian*[ones(Np,1);zeros(Nr,1)])
    
    eigenvalues=study_eigenvalues(out, '(prec)^{-1}J',1);
    %eig6=eigs(@(x) jacobian*x,Nr+Np,40)
    
    figure
    plot(eigenvalues,'o')

    return
  end
  
  outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) apply_saddle_point(x,@(y) A*y ,...
																														@(y) B1T*proj(y),...
																														@(z) B2*z,...
																														@(z) C*z,...
																														Np,Nr), ...
																		rhs, ctrl_outer, @(z) prec(z),[],controls.left_right );
	for i = 1:N
		imb=d(Np+1+(i-1)*ncellrho:Np+i*ncellrho)'*area2h;
		state_message=sprintf('y*area= %1.4e',imb);
		fprintf('%s \n',state_message);
	end

	for i = 1:N+1
		imb=sum(F.p(1+(i-1)*ncellphi:i*ncellphi));
		state_message=sprintf('imb= %1.4e',imb);
		fprintf('%s \n',state_message);
	end

	
  if controls.remove_imbalance
    % restore solution
    increment=zeros(Np+Nr+N,1);
    increment(1:Np)        = d(1:Np);
    increment(Np+1:Np+Nr+N)= mat_H*d(Np+1:Np+Nr);
    d=increment;
    Nr=Nr+N;

  end
  
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


  % assembly full system
  dp = d(1:Np); dr = d(Np+1:end);
  ds = JF.ss\(-F.s-JF.sr*dr);
  d = [dp; dr; ds];
  [res,resp,resr,ress]=compute_linear_system_residuum(JF,F,d)

  relres=res;
 

  
  normd = norm(d);

  if (0)    
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
      
      
      %saveas(gcf,strcat(controls.basename,controls.sys_name,'denseS.png'));
      saveas(gcf,'denseS.png');
      
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

  if (1)
    out=prec_times_matrix(prec,jacobian);
    
    eigenvalues=study_eigenvalues(out, '(prec)^{-1}J',1);
    %eig6=eigs(@(x) jacobian*x,Nr+Np,40)
    
    figure
    plot(eigenvalues,'o')

    return
  end


  
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

  if (1)    
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

  if (strcmp(controls.gamma,'auto'))
    gamma=norm(A,'fro')/(norm(B1T,'fro')^2);
  else
    gamma=controls.gamma;
  end				% init diag(C)^{-1}

  
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
  elseif ( strcmp(W_approach,'Ones') )
    W=speye(Nr);
    invW=speye(Nr);%sparse(1:Nr,1:Nr,(1.0./spdiags(W,0))',Nr,Nr);
  elseif ( strcmp(W_approach,'select') )
    W=speye(Nr);
    selection=zeros(Nr,1);
    for i=1:N
      selection(1+(i-1)*ncellrho)=1.0;
    end
    invW=sparse(1:Nr,1:Nr,(selection)',Nr,Nr);
  elseif ( strcmp(W_approach,'rho') )
    W=speye(Nr);
    invW=-JF.ss;
  elseif ( strcmp(W_approach,'rank') )
    W=speye(Nr);
    vectors=zeros(Nr,N);
    for i=1:N
      vectors(1+(i-1)*ncellrho:i*ncellrho,i)=1;
    end
    invW=vectors*vectors';
  end
    

  build_inv_augA=tic;
				% augemented system

  %figure
  %plot(selection)
  augA=A+gamma*B1T*invW*B2;
  %figure
  %spy(augA)

  %eigenvalues=study_eigenvalues(full(augA),'augA',1);%controls.logID);
  %eigenvalues(Np-8:Np-1)
  %eigenvalues(Np)
  %eigenvalues(1)

  %figure
  %plot(eigenvalues,'o')

  %return
  augB1T=B1T*(speye(Nr,Nr)-gamma*invW*C);
  %norm(nonzeros(augB1T))
  aug_f=rhs(1:Np)+gamma*B1T*invW*rhs(Np+1:Np+Nr);
  augrhs=zeros(Np+Nr,1);
  augrhs(1:Np)=aug_f; %f_gamma
  augrhs(1+Np:Np+Nr)=rhs(Np+1:Np+Nr); %g

  % init inverse of augA
  
  tA=tic;

  approach_inverse_A=controls.approach_inverse_A;
  if ( strcmp(approach_inverse_A,'full'))
    ctrl_loc=ctrl_solver;
    index_agmg=0;
    index_agmg=index_agmg+1;
    ctrl_loc.init(controls.ctrl_inner11.approach,...
		  controls.ctrl_inner11.tolerance,...
		  controls.ctrl_inner11.itermax,...
		  controls.ctrl_inner11.omega,...
		  controls.ctrl_inner11.verbose,...
		  'augA',index_agmg);

    inv_augA=sparse_inverse;
    inv_augA.init(augA,ctrl_loc);
    inv_augA.cumulative_iter=0;
    inv_augA.cumulative_cpu=0;
    inv_augA.dimblock=ncellphi;
				%inv_S.info();
    
    inverse11 = @(x) inv_augA.apply(x);

    inner_nequ=Np;
    
  elseif( strcmp(approach_inverse_A,'block'))
    index_agmg=0;
    index_agmg=index_agmg+1;
    
    % define array of sparse inverse
    diag_block_invaugA(Nt,1)=sparse_inverse;

    % create blocks
    nAi=ncellphi;
    ctrl_loc=ctrl_solver;
    for i=1:Nt
      index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(controls.ctrl_inner11.approach,...
		    controls.ctrl_inner11.tolerance,...
		    controls.ctrl_inner11.itermax,...
		    controls.ctrl_inner11.omega,...
		    controls.ctrl_inner11.verbose,...
		    sprintf('%s Aaug%d',controls.ctrl_inner11.label,i),...
		    index_agmg);
      diag_block_invA(i).name=sprintf('inverse A%d',i);
      
      % get block add relaxation
      matrixaugAi=augA((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      
      % define inverse
      diag_block_invaugA(i).init( matrixaugAi,ctrl_loc);     

    end
    % define function
    inverse11  = @(x) apply_block_diag_inverse(diag_block_invaugA,x);
  end
  

  
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
  elseif( strcmp(S_approach,'Adiag'))
			       % define explicitely inverse of diag(A)
    invDiagA = sparse(1:Np,1:Np,(1.0./spdiags(augA,0))',Np,Np);

    
    % assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
    if (    strcmp(W_approach,'rank'))
      approx_SCA = sparse(C+B2*invDiagA*B1T);
    else
      approx_SCA = sparse(C+B2*invDiagA*augB1T);
    end
    
    % assembly approximate inverse
    inv_S=sparse_inverse;


    % invDiagA)we need a new seed for agmg, in case is used
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
    
    inv_S.init(approx_SCA+controls.relax4inv22*speye(Nr,Nr),ctrl_loc);
    inv_S.info_inverse.label='schur_ca';
    inv_S.cumulative_iter=0;
    inv_S.cumulative_cpu=0;
    
    inverse_cpu=inv_S.init_cpu;
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

  
  if (0)
    % prec_ex=[augA augB1T; sparse(Nr,Np) -approx_SCA];
    % inv_P=sparse_inverse;
    % ctrl_loc.init('direct',...
    % 		  ctrl_inner22.tolerance,...
    % 		  ctrl_inner22.itermax,...
    % 		  ctrl_inner22.omega,...
    % 		  ctrl_inner22.verbose,...
    % 		  'P^-1',...
    % 		  index_agmg);
    
    % inv_P.init(prec_ex,ctrl_loc);
    % prec_ex=@(x) inv_P.apply(x);
    %eigenvalues=study_eigenvalues(full(augjacobian), '(prec)^{-1}J',1);
    %eig6=eigs(@(x) jacobian*x,Nr+Np,40)
    
    %figure
    %plot(eigenvalues,'o')

    

    out=prec_times_matrix(prec,augjacobian);

    
    
    eigenvalues=study_eigenvalues(out, '(prec)^{-1}J',1);
    %eig6=eigs(@(x) jacobian*x,Nr+Np,40)
    
    figure
    plot(eigenvalues,'o')

    return
  end

  
  if ( verbose >= 2)
    fprintf('START SOLVER \n')
  end
  [d,info_J]=apply_iterative_solver(@(x) augjacobian*x, augrhs, ctrl_outer, prec );
  if ( verbose >= 2)
    fprintf('END SOLVER \n')
  end
  %info_J.print();
  outer_iter=uint64(info_J.iter);
  outer_cpu=toc(outer_timing);

  if ( strcmp(approach_inverse_A,'full'))
    inner_iter=inv_augA.cumulative_iter;
  elseif( strcmp(approach_inverse_A,'block'))
    inner_iter=0;
    for i=1:Nt
      inner_iter=inner_iter+diag_block_invaugA(i).cumulative_iter;
    end
      
  end
  
  
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
  if ( strcmp(approach_inverse_A,'full'))
    inv_augA.kill();
  elseif( strcmp(approach_inverse_A,'block'))
    for i=1:Nt
      diag_block_invaugA(i).kill();
    end
  end
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
  
elseif (sol == 17)
  % get variable for assembling system
  rho=-spdiags(JF.ss);
  deltat=1/(N+1);

	% assembly lambda residuum
	% f(k)=(area*rho_k-1)
	v = JF.area2h;
	f_lambda = zeros(N,1);
	W_mat=zeros(N,Nr);
	for k = 1:N
		f_lambda(i) = v' * get_slice(rho,k, ncellrho, N) - 1; % F^k_lambda
    W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = v'; % d F^k_lambda drho_k 
  end
  
  % build projector matrix
	factor=1/(norm(v))^2
  P=speye(Nr,Nr)-factor*W_mat'*W_mat;
  rhs = -[F.p; F.r; F.s; f_lambda];

  % hessian + sostituion for slack varaible
  Sys = sparse([...
		 JF.pp,         JF.pr,  JF.ps,        sparse(Np,N) ;
		 JF.rp,         JF.rr,  JF.rs,        W_mat'        ;
		 sparse(Nr,Np), JF.sr,  JF.ss,        sparse(Nr,N);
		 sparse(N, Np), W_mat,  sparse(N,Nr), sparse(N, N) ...
	       ]);

	
	
  % solve
  dim=Np+Nr+Nr+N;
	%Sys(indc,:) = sparse(1,dim); Sys(:,indc) = sparse(dim,1); Sys(indc,indc) = 1; rhs(indc)=0; 
	increment=Sys\rhs;
  
  % get phi,rho,s increment and lambda
  d=increment(1:Np+Nr+Nr);
  dl=increment(1+Np+Nr+Nr: Np+Nr+Nr+N);
	[resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
	resume_msg=sprintf('res=%1.1e | p=%1.1e r=%1.1e s=%1.1e\n',resnorm,resp,resr,ress);
	disp(resume_msg)

	compare=0;
	if (compare)
		indc=1
		SysR = sparse([...
  									JF.pp,         JF.pr,   JF.ps, ;
  									JF.rp,         JF.rr,   JF.rs;
  									sparse(Nr,Np), JF.sr,  JF.ss   ...
  								]);
		rhs = -[F.p; F.r; F.s];
		dim=Np+2*Nr;
		%SysR(indc,:) = sparse(1,dim); SysR(:,indc) = sparse(dim,1); SysR(indc,indc) = 1; rhs(indc)=0; 
		exact_d=SysR\rhs;
		error=d-exact_d;
		figure
		plot(error(1:Np+Nr))
	end

	
	d(1:Np) = dp_correction(d(1:Np),JF,dl);

	if (compare)
		error=d-exact_d;
		hold on
		plot(error(1:Np+Nr))
		hold off
	end
  [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
	resume_msg=sprintf('res=%1.1e p=%1.1e r=%1.1e s=%1.1e\n',resnorm,resp,resr,ress);
	
  normd = norm(d);
 
 
    
  inner_nequ=0;
  inner_iter=0;
  outer_iter=0;
  outer_cpu=0;
  prec_cpu=0;
  inner_nequ2;
  inner_iter2=0;
  inner_nequ3=0;
  inner_iter3=0;


elseif (sol==18)
  % get variable for assembling system
  rho=-spdiags(JF.ss);
  deltat=1/(N+1);
  
  % build matrix W
  W_mat=zeros(N,Nr);
	v=deltat*area2h';
  for i = 1:N
    W_mat(i,1+(i-1)*ncellrho:i*ncellrho)=v';
  end

	
	factor=1/(norm(v))^2
  P=speye(Nr,Nr)-factor*W_mat'*W_mat;

  % assembly lambda residuum
  
  f_lambda=zeros(N,1);
  for i = 1:N
    f_lambda(i)=deltat*(area2h' * rho(1+(i-1)*ncellrho:i*ncellrho)-1);
		fprintf('flambda(%d)=%1.2e\n',i,f_lambda(i)) 
  end

	SysR = sparse([...
  								JF.pp,         JF.pr,   JF.ps, ;
  								JF.rp,         JF.rr,   JF.rs;
  								sparse(Nr,Np), JF.sr,  JF.ss   ...
  							]);
	rhs = -[F.p; F.r; F.s];

	exact_d=SysR\rhs;
	res=norm(SysR*exact_d - rhs);
	fprintf('res real=%1.2e\n',res)

	print_imbalance(F.p,ncellphi)
	norm(P'*exact_d(Np+1:Np+Nr)-exact_d(Np+1:Np+Nr))
	disp('here')
	for i = 1:N
		imb=exact_d(Np+1+(i-1)*ncellrho:Np+i*ncellrho)'*area2h;
		state_message=sprintf('y*area= %1.4e',imb);
		fprintf('%s \n',state_message);
	end
	disp('here')
	W_mat*exact_d(Np+1:Np+Nr)

	disp('proj')
	if (1)
		rhsR = -[F.p; P*F.r; F.s];
		SysR = sparse([...
  									JF.pp,         JF.pr,   JF.ps, ;
  									P*JF.rp,       P*JF.rr, P*JF.rs;
  									sparse(Nr,Np), JF.sr,  JF.ss   ...
  								]);
		
		
		d = SysR\rhsR;
		res=norm(SysR*d - rhsR);
		fprintf('res aug=%1.2e\n',res)
	else
		rhsR = [f1; P*f2];
		SysR = sparse([...
										A,   B1T*P;
										P*B2, -P*C*P;
									]);
		
		d = SysR\rhsR;
		res=norm(SysR*d - rhsR);
		ds = Dr\(h-Ds*d(Np+1:Np+Nr));
		d = [d; ds];
		fprintf('res aug=%1.2e\n',res)
	end

	% get lambda increment dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)
  %dl=get_dl(JF,F,d);
  
  % correction of phi increment 
	%d(1:Np) = dp_correction(d(1:Np),JF,dl);

	%figure
	%error=d(1:Np+Nr)-exact_d(1:Np+Nr);
	%plot(error)
  
  % dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z-g)
  res_r=-([JF.rp JF.rr JF.rs]*d + F.r);

	disp(norm(res_r))
  factor=1/((deltat*norm(area2h))^2);
  dl=factor* W_mat*(res_r);
  %dl=get_dl(JF,F,d)
  for i = 1:N    
    fprintf('dlambda(%d)=%1.2e\n',i,dl(i)) 
  end

	rhsP = -[F.p; F.r; F.s; f_lambda];

  % hessian + sostituion for slack varaible
  SysP = sparse([...
		 JF.pp,         JF.pr,  JF.ps,        sparse(Np,N) ;
		 JF.rp,         JF.rr,  JF.rs,        W_mat'        ;
		 sparse(Nr,Np), JF.sr,  JF.ss,        sparse(Nr,N);
		 sparse(N, Np), W_mat,  sparse(N,Nr), sparse(N, N) ...
								]);

	resP=SysP*[d;dl]-rhsP;
	norm(SysP*[d;dl]-rhsP)
	figure
	plot(resP)

	return
  
  %correction of phi increment 
  sum_dl=0;
  for i = 1:N
    sum_dl=sum_dl+deltat*dl(i);
		j=i;
    d(1+(j-1)*ncellphi:j*ncellphi)=d(1+(j-1)*ncellphi:j*ncellphi)+sum_dl;
  end

	

	hold on
	error=d(1:Np+Nr)-exact_d(1:Np+Nr);
	plot(error)
  hold off
  
  [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);

	resume_msg=sprintf('res=%1.1e p=%1.1e r=%1.1e s=%1.1e\n',resnorm,resp,resr,ress);
	disp(resume_msg)
	
  normd = norm(d);
 
	return
    
  inner_nequ=0;
  inner_iter=0;
  outer_iter=0;
  outer_cpu=0;
  prec_cpu=0;
  inner_nequ2;
  inner_iter2=0;
  inner_nequ3=0;
  inner_iter3=0;

elseif (sol == 19)
  % get variable for assembling system
  rho=-spdiags(JF.ss);
  deltat=1/(N+1);
  
  % build matrix W
  W_mat=zeros(N,Nr);
  for i = 1:N
    W_mat(i,1+(i-1)*ncellrho:i*ncellrho)=deltat*area2h';
  end

  %P=speye(Nr,Nr)-W_mat'*W_mat;

	% define implicetly projector
	factor=1/(deltat*norm(area2h))^2 
	proj=@(rho) rho - factor * W_mat'*(W_mat*rho);
  
  % assembly lambda residuum
  
  f_lambda=zeros(N,1);
  for i = 1:N
    f_lambda(i)=deltat*(area2h' * rho(1+(i-1)*ncellrho:i*ncellrho)-1);
		fprintf('flambda(%d)=%1.2e\n',i,f_lambda(i)) 

  end
  
  rhsR = [f1; proj(f2)];
 
	%
	% PREPROCESS TIME
	%
  time_prec=tic;
	time_A=tic;
	
  % define approxiamte inverse of (~A)^{-1}  
  inverseA_approach=controls.inverse11;
	ctrl_inner11 = controls.ctrl_inner11;

	disp(inverseA_approach)
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    invA=sparse_inverse;
    invA.init(A+controls.relax4inv11*speye(Np,Np),controls.ctrl_inner11);
    invA.cumulative_iter=0;
    invA.cumulative_cpu=0;

    % define function
    inv_A = @(x) invA.apply(x);
  elseif( strcmp(inverseA_approach,'diag'))
    % partion matrix 
    nAi=ncellphi;
   
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
    		    sprintf('BD%sA%d',ctrl_inner11.label,i));
      diag_block_invA(i).name=sprintf('inverse A%d',i);

      % create local block and passing to solver
      % with a potential relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + ...
	       controls.relax4inv11*speye(nAi,nAi) ;
      diag_block_invA(i).init(matrixAi,ctrl_loc);     
    end

    % define function
    inv_A  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end
	build_inverse11=toc(time_A);
	

  %
  % BUILD APPROXIMATE INVERSE OF SCHUR COMPLEMENT
  %
  
	% set dimension of inner solver
  inner_nequ=Nr;

	% store time
  build_S=tic;

	%select among differnt approaches
	inverse22 = controls.inverse22;
	ctrl_inner22 = controls.ctrl_inner22;
	
  if (strcmp(inverse22,'diagA'))
    % assembly SAC=-(C+B2 * diag(A)^{-1} B1T)
    invDiagA   = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);
    approx_SCA = (C+B2*invDiagA*B1T);

    inv_SCA=sparse_inverse;
    inv_SCA.init(approx_SCA + controls.relax4inv22*speye(Nr,Nr),ctrl_inner22);
    inv_SCA.info_inverse.label='schur_ca';
    inv_SCA.cumulative_iter=0;
    inv_SCA.cumulative_cpu=0; 
    inverse_block22 = @(x) -inv_SCA.apply(x);

  elseif (strcmp(inverse22,'diagA'))
    % define the preconditioner as the result of an iterative solver
    % acting with operator C+B^T A^{-1} B

    SCA = @(x) (C*x + B2*(invA(B1T*x)));
    %assembly inverse of  SAC iterative solver with no precodnitioner
    inverse_block22 = @(y) ...
											- apply_iterative_solver( SCA, y, ctrl_inner22, @(z) z);
  end
  cpu_assembly_inverseS=toc(build_S);

  % store time required to build prec
  prec_cpu=toc(time_prec);

  %
  % Define action of preconditoner
  % 
  prec = @(x) SchurCA_based_preconditioner(x, ...
																					 @(y) projector(inv_A(y),kernel),...
																					 @(z) proj(inverse_block22(z)),...
																					 @(y) B1T*proj(y),...
																					 @(z) proj(B2*z),...
																					 Np,Nr,...
																					 controls.outer_prec,ncellphi);

	outer_timing=tic;
  [d,info_J]=apply_iterative_solver(@(x) ...
																		 apply_saddle_point(x,@(y) A*y ,...
																												@(y) B1T*proj(y),...
																												@(z) proj(B2*z),...
																												@(z) proj(C*proj(z)),...
																												Np,Nr), ...
																		rhsR, controls.ctrl_outer, ...
																		@(z) prec(z),[],controls.left_right );

	info_J.print();

	
	print_imbalance(d(1:Np),ncellphi)

	
	for i = 1:N
		imb=d(Np+1+(i-1)*ncellrho:Np+i*ncellrho)'*area2h;
		state_message=sprintf('y*area= %1.4e',imb);
		fprintf('%s \n',state_message);
	end

	res=apply_saddle_point(d,@(y) A*y ,...
														@(y) B1T*y,...
														@(z) B2*z,...
														@(z) C*z,...
														Np,Nr)-rhsR;
	
	fprintf('|res p|=%1.1e - |res r|=%1.1e - |rhsR|=%1.1e  \n',...
					norm(res(1:Np)),norm(res(1+Np:Np+Nr)),norm(rhsR));

	d(Np+1:Np+Nr)=proj(d(Np+1:Np+Nr));
	res=apply_saddle_point(d,@(y) A*y ,...
												 @(y) B1T*y,...
												 @(z) B2*z,...
												 @(z) C*z,...
												 Np,Nr)-rhsR;
	
	fprintf('|res p|=%1.1e - |res r|=%1.1e - |rhsR|=%1.1e  \n',...
					norm(res(1:Np)),norm(res(1+Np:Np+Nr)),norm(rhsR));
	
	
	% get s increment
  ds = JF.ss\(-F.s-JF.sr*d(Np+1:end));
  d = [d; ds];
	

	
	
	% stop cpu
  outer_cpu=toc(outer_timing);

	% get info
  inner_iter=inv_SCA.cumulative_iter;
  outer_iter=uint32(info_J.iter);
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;
  normd = norm(d);

	% get lambda increment dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)
  dl=get_dl(JF,F,d);
  
  % correction of phi increment 
	d(1:Np) = dp_correction(d(1:Np),JF,dl);
  [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
  
  resume_msg=sprintf('res=%1.1e p=%1.1e r=%1.1e s=%1.1e\n',resnorm,resp,resr,ress);
  normd = norm(d);

	disp(resume_msg)
	

	return
 
    
  inner_nequ=Np;
  inner_iter=0;
  outer_iter=0;
  outer_cpu=0;
  prec_cpu=0;
  inner_iter2=Nr;
  inner_nequ3=0;
  inner_iter3=0;

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
