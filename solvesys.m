function [d, sol_stat, normd,resume_msg] = solvesys(JF,F,controls)

% solve the linear system d = -JF\F, reducing the system by eliminating the
% additional variable s and imposing the condition d(indc) = 0 to fix the 
% constant of the potential dp
  
sol=controls.sol;

resume_msg='';

% counter for agmg solver.
% see inverse_sparse why we need this
index_agmg=0;

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
swap_sign=1;

if (swap_sign)
  sign=-1;
else
  sign=1;
end

				% reduced system to phi rho variables
A = sparse(JF.pp);
B1T = sparse(JF.pr-JF.ps*(JF.ss\JF.sr));
B1T_time  = JF.B1T_time;
B1T_space = JF.B1T_space;
B2 = sparse(JF.rp);
R = sparse(JF.rr);
M   = sparse(JF.rs);
Ds = sparse(JF.sr);
Dr = sparse(JF.ss);
C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
temp=(JF.ss\h);
f1 = f-JF.ps*temp;
f2 = g-JF.rs*temp;

Mphi=JF.Mxt;
Mrho=-JF.rs;
inv_Mphi  = sparse(1:Np,1:Np,(1./spdiags(Mphi,0))',Np,Np);
inv_Mrho  = sparse(1:Nr,1:Nr,(1./spdiags(Mrho,0))',Nr,Nr);

% swap C sign for having standard saddle point notation 
C = -C;
if (swap_sign)
  A=-A;
  B1T=-B1T;
	B1T_time =-B1T_time;
	B1T_space=-B1T_space;
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

% set the vectors in the kernels of A
A_kernels=zeros(N+1,Np);
for k=1:N+1
	 A_kernels(k,1+(k-1)*ncellphi:k*ncellphi)=ones(1,ncellphi)/sqrt(ncellphi);
end

area2h=JF.area2h;

deltat=1/(N+1);
  
if sol==1
  % Solve using Matlab's backslash  
	if (controls.ground == 1)
		% Ground matrix A, B1T and f1.
		% Change the upper bound cycle to see the effect of
		% "grounding" all blocks. 
		for iblock = 1:1
			irow = controls.ground_node+(iblock-1)*ncellphi;
			A(irow,:) = sparse(Np,1);
			A(irow,irow)= 1.0;
			B1T(irow,:) = sparse(Nr,1);
		end
	else
		% Other relaxation type
		A=A+contols.relax11*speye(Np,Np);
	end

	% we can solve the full, the partially reduce or the fully reduced system
	outer_cpu=tic;
	if (controls.reduced == 0)
		% form the full matrix and solve it
		Sys = sparse([A B1T sparse(Np,Nr); B2 R M; sparse(Nr,Nr) Ds Dr]);
		rhs = [f;g;h]
		d = Sys\rhs;

	elseif (controls.reduced == 1)
		% form the saddle point system removing the slackness variable
		Sys = sparse([A B1T; B2 -C ]);
		rhs = [f;f2];
		d = Sys\rhs;
		% get the ds 
		ds = Dr\(h-Ds*d(Np+1:Np+Nr));
		d=[d;ds];

		outer_iter=0;
		inner_iter=0;
		
		
	elseif (controls.reduced == 2)
		% form the primal Schur Complement, solve w.r.t to dp
		inv_C = sparse(1:Nr,1:Nr,(1.0./spdiags(C,0))',Nr,Nr);
		Sys = A+B1T*inv_C*B2;
		rhs = f1-B1*(inv_C*f2); fp = fp(2:end);

		% get dr and ds
    dr = inv_C*(-f2-B2*dp);
    ds = Dr\(h-Ds*dr);
    d = [dp; dr; ds];
	end
	prec_cpu=0;
	normd=norm(d);
	outer_cpu=toc(outer_cpu);


% elseif sol==100
	
		
	
%   if (0)
%     Sys = [A B1T; B2 -C];
%     fp=f1;

% 		indc=1;
% 		%[dmax,indc]=max(spdiags(A(1:ncellphi:1:ncellphi),0));
%     Sys(indc,:) = sparse(1,Np+Nr);
% 		Sys(:,indc) = sparse(Np+Nr,1);
% 		Sys(indc,indc) = 1;
% 		fp(indc)=0; 
%     d = Sys\[fp;f2];

		

% 		[A,B1T,B2,C,f1,f2,W_mat]=get_saddle_point(JF,F);
% 		B1T_time=-JF.B1T_time;
% 		B1T_space=-JF.B1T_space;

% 		P = @(rho) ort_proj(rho,W_mat);
% 		B1T_time=matrix_times_prec(B1T_time,P);
% 		B1T_space=matrix_times_prec(B1T_space,P);
% 		%C=matrix_times_prec(C,P)+1e-12*speye(Nr,Nr);
% 		[A,B1T_time,B1T_space,B2,f1]=ground_saddle(A,B1T_time,B1T_space,B2,f1,N);
% 		B1T=B1T_time+B1T_space;


		
% 		diagonal_scaling=0;
% 		if (diagonal_scaling)
% 			% scale system by its diagonal
% 			[sqrt_diagA,sqrt_diagC,inv_sqrt_diagA,inv_sqrt_diagC] = get_diagonal_scaling (A,C);
			
% 			[A,B1T,B2,C,fp,f2,W_mat]=scale_saddle_point(A,B1T,B2,C,...
% 																									fp,f2,W_mat,...
% 																									inv_sqrt_diagA,...
% 																									inv_sqrt_diagC);
% 		end
% 		Sys = [A B1T; B2 -C];
% 		d_new=Sys\[fp;f2];
		
		

% 		if (diagonal_scaling)
% 			d_new(1:Np)=inv_sqrt_diagA*d_new(1:Np);
% 			d_new(1+Np:Np+Nr)=inv_sqrt_diagC*d_new(1+Np:Np+Nr);
% 		end
% 		%d_new(1+Np:Np+Nr)=P(d_new(1+Np:Np+Nr));
	

	

% 		 E=zeros(Np,N);
% 		 for k=1:N
% 		 	E(k*ncellphi+1:(k+1)*ncellphi,k)=1.0;
% 		 end
	

% 		 W_mat=zeros(N,Nr);
% 		 for k = 1:N
% 		 	 W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = JF.area2h';
% 		 end

% 		 MM=B2*E;
		 
% 		 matrix=C\(MM);
% 		 %figure;
% 		 %spy(matrix)
% 		 matrix_r=W_mat*(matrix)
% 		 rhs_r=W_mat*d_new(Np+1:Np+Nr)
% 		 W_mat*d(Np+1:Np+Nr)
		
% 		 alphas=matrix_r\rhs_r

% 		 figure
% 		 plot(d_new(Np+1:Np+Nr)-d(Np+1:Np+Nr))

% 		 dp=d_new(1:Np)+E*alphas;
% 		 dr=d_new(Np+1:Np+Nr)+matrix*alphas;

		 

% 		 ds = Dr\(h-Ds*d_new(Np+1:Np+Nr));
%      d_new = [d_new;ds];
% 		 [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d_new)

% 		 return

		 
	
%     res=norm(Sys*d-[f1;f2])/norm([f1;f2]);
%     dp = d(1:Np); dr = d(Np+1:end);
%     ds = Dr\(h-Ds*dr);
%     d = [dp; dr; ds];

%     normd = norm(d);
%     resume_msg=sprintf('direct solver res=%1.2e',norm(res));
    
    
%   else
%     Sys = sparse(...
% 		  [JF.pp,        JF.pr, JF.ps;
% 		   JF.rp,        JF.rr, JF.rs;
% 		   sparse(Nr,Np), JF.sr, JF.ss]...
% 	    );
    
%     rhs=-[f;g;h];
%     Sys(indc,:) = sparse(1,Np+Nr+Nr); Sys(:,indc) = sparse(Np+Nr+Nr,1); Sys(indc,indc) = 1; rhs(indc)=0.0;
%     d = Sys\rhs;
%     res=norm(Sys*d-rhs)/norm(rhs);


    
%     %eigenvalues= study_eigenvalues(full(Sys),'J_Full');
%     %figure
%     %plot(eigenvalues,'o')
%     %return

% 				%dp = d(1:Np); dr = d(Np+1:end);
%     %ds = Dr\(h-Ds*dr);
%     %d = [dp; dr; ds];

%   end

%   if (0)
%     % we enforse the rho increment to have zero mass

    
%     % P=  (I_{m-1}
%     %    -(a1 -a2 -a_{m-1})/a_m 
%     mat_P=speye(ncellrho,ncellrho-1);   
%     mat_P(ncellrho,1:ncellrho-1)=-area2h(1:ncellrho-1)'/area2h(ncellrho);

%     % H=block diagonal (P) 
%     mat_H = repmat({mat_P},1,N);
%     mat_H = blkdiag(mat_H{:});


%     s_B1T = B1T * mat_H;
%     s_B2 = mat_H' * B2;
%     s_C  = mat_H' * C * mat_H;

%     s_rhs = [f1; mat_H' * f2];

    
%     J_short=[A, s_B1T; s_B2, -s_C];

%     d_short = J_short\s_rhs;

%     dp = d_short(1:Np);
%     dr = mat_H*d_short(Np+1:size(d_short,1));    
%     ds = Dr\(h-Ds*dr);



%     res_r=-(JF.rp * dp + JF.rr * dr + JF.rs * ds + F.r);
%     dl=zeros(N,1);

%     deltat=1/(N+1);
%     factor=1/(deltat*norm(area2h))^2;
%     for k = 1:N    
%       dl(k)= factor * (area2h'* (res_r(1+(k-1)*ncellrho:k*ncellrho))) ;
%       fprintf('dlambda(%d)=%1.2e\n',k,dl(k)) 
%     end
    
%     sum_dl=0;

% 		disp(i)
%     for k = 1:N
%       sum_dl=sum_dl+dl(k);
%       dp(1+(k-1)*ncellphi:k*ncellphi)=dp(1+(k-1)*ncellphi:k*ncellphi)+deltat*sum_dl;
%     end
    
%     d = [dp; dr; ds];

%     [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
      
    
% 				% redefine rho dimension 

%   end


%   [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
%   normd = norm(d);
%   resume_msg = sprintf('direct solver res=%1.2e',resnorm);

%   %figure
%   %spy(Sys)
%   %return
    
%     inner_nequ=0;
%     inner_iter=0;
%     outer_iter=0;
%     outer_cpu=0;
%     prec_cpu=0;
%     inner_nequ2;
%     inner_iter2=0;
%     inner_nequ3=0;
%     inner_iter3=0;
    
elseif sol==2
    
    % Solve with the Schur complement A|C, using Matlab's backslash (or agmg)
  bs=tic;
    A = JF.pp; B1 = JF.pr; B2 = JF.rp; inv_C = -((JF.sr*JF.rs)\JF.ss);
    f1 = F.p;
    f2 = F.r-JF.rs*(JF.ss\F.s);
    S = -A+B1*(inv_C*B2); S = S(2:end,2:end);
    fp = f1-B1*(inv_C*f2); fp = fp(2:end);
    dp = S\fp;
    outer_cpu=toc(bs);
    outer_iter=0;
    sol_stat = struct('flag',0,'relres',[],'iter',[]);
    
    %[dp,flag,res,iter] = agmg(S,fp,0,1e-6,100,0,[],0);
    %sol_stat = struct('flag',flag,'relres',res,'iter',iter);

    inner_iter=0;
    dp = [0; dp];
    dr = inv_C*(-f2-B2*dp);
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
  inv_C = sparse(1:Nr,1:Nr,(1.0./spdiags(C))',Nr,Nr);
 
  % assembly schur AC S=A+B1
  % reduced to system only in phi
  S = A+B1T*(inv_C*B2);
  fp = f1+B1T*(inv_C*f2);

  %norm(S*d(1:Np)-fp)
  
  % ground the solution
  S(indc,:) = sparse(1,Np);
  S(:,indc) = sparse(Np,1);
  S(indc,indc) = 1;
  fp(indc)     = 0;
  
    
  % solve with multigrid
  %[dp,flag,res,iter] = agmg(S,fp,0,ctrl.tolerance,ctrl.itermax,0,[],0);
  inv_S=sparse_inverse;
  inv_S.init(S,ctrl_inner11);
  dp=inv_S.apply(fp);
  %inv_S.info.print();
  flag=inv_S.info.flag;
  relres=inv_S.info.res;
  iter=inv_S.info.iter;
  %inv_S.kill();

  %res;norm(S*dp-fp)

  % get solution of system JF d  = F  
  dr = inv_C*(B2*dp-f2);

  %J=[A, B1T;  B2, -C];
  %norm(J*[dp;dr]-[f1;f2])

  ds = JF.ss\(h-JF.sr*dr);
  d = [dp; dr; ds];
    
  normd = norm(d);

% elseif sol==7
%   % rhs
%   f=-F.p;
%   g=-F.r;
%   h=-F.s;
  
%   % Solve wsith the Schur complement A|C, using Matlab's backslash (or agmg)
%   swap_sign=1;

%   % reduced system to phi rho variables
%   A = JF.pp; B1T = JF.pr; B2 = JF.rp;
%   C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
%   f1 = f;
%   f2 = g-JF.rs*(JF.ss\h);
  
%   % swap C sign for having standard saddle point notation 
%   C = -C;
%   if (swap_sign)
%     A=-A;
%     B1T=-B1T;
%     B2=-B2;
%     C=-C;
%     f1=-f1;
%     f2=-f2;
%   end

%   %S=compute_SchurCA(A,B1T,B2,C,Nt,ncellrho,ncellphi)
%   %S=A+B1T*(C\B2);
  
%   %J=[A, B1T;  B2, -C];
%   %norm(J*d(1:Np+Nr)-[f1;f2])

%   % start solution
%   % invert C
%   diagC=spdiags(C);
%   %fprintf('(%9.4e <= C <= %18.12e) \n',min(diagC),max(diagC))
%   inv_diagC = sparse(1:Nr,1:Nr,(1.0./diagC)',Nr,Nr);
 
%   % grounding
%   A(indc,:)= sparse(1,Np); %A(:,indc)= sparse(Np,1);
%   A(indc,indc)=1;
%   B1T(indc,:)= sparse(1,Nr);
%   rhs(indc)=0;

  
%   % assembly schur AC S=A+B1
%   % reduced to system only in phi
%   S = A+B1T*(inv_diagC*B2);
%   fp = rhs(1:Np)+B1T*(inv_diagC*rhs(1+Np:Np+Nr));

%   if ( compute_eigen )
%     disp("START EIG")
%     eigenvalues=eig(full(S));
%     disp("END EIG")
%     for i=1:4
%       fprintf('(%9.4e <= eigvalues <= %18.12e) \n',eigenvalues(i),eigenvalues(Np-i+1))
%     end
%   end
    
%   % solve

%   % init solver
%   inv_S=sparse_inverse;
%   inv_S.init(S,ctrl_inner11);
% 				% apply inverse
%   dp=inv_S.apply(fp);
%   inner_nequ=inv_S.nequ;
  
%   % copy info
%   flag=inv_S.info_inverse.flag;
%   relres=inv_S.info_inverse.res;
%   inner_iter=inv_S.info_inverse.iter;
%   outer_iter=0;
%   outer_cpu=0.0;
%   %inv_S.kill();
  
%   % get solution of system JF d  = F  
%   dr = invdiagC*(B2*dp-rhs(1+Np:Np+Nr));

%   %J=[A, B1T;  B2, -C];
%   %norm(J*[dp;dr]-[f1;f2])

%   ds = JF.ss\(h-JF.sr*dr);
%   d = [dp; dr; ds];
    
%   normd = norm(d);

% elseif sol==8
%   % Invert with respect to C block and try to "invert" S=A+BT*C^{-1} B
 
%   % start solution
%   % invert C
%   fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)))
%   inv_C = sparse(1:Nr,1:Nr,(1.0./spdiags(C))',Nr,Nr);

%   % grounding
%   A(indc,:)= sparse(1,Np);
%   A(indc,indc)=1;
%   B1T(indc,:)= sparse(1,Nr);
%   rhs(indc)=0;

  
%   % assembly schur AC S=A+B1
%   % reduced to system only in phi
%   S = A+B1T*(inv_C*B2);
%   fp = rhs(1:Np)+B1T*(inv_C*rhs(1+Np:Np+Nr));

%   %
%   % partion diagonal and upper diagonal blocks
%   %
%   diag_block_S       = cell(Nt, 1);
%   below_diag_block_S= cell(Nt-1, 1);
%   nAi=ncellphi;
%   for i=1:Nt-1
%     diag_block_S{i}      =S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
%     below_diag_block_S{i}=S(i    *nAi+1 : (i+1)*nAi , (i-1)*nAi+1 : i*nAi);
%   end
%   i=Nt;
%   diag_block_S{i}=S((i-1)*nAi+1:i*nAi,(i-1)*nAi+1:i*nAi);

%   %				
%   % create define inverse on diagonal block
%   % relax diagonal by a factor 1/omega 0<omega<1
%   %
%   diag_block_inv_S(1:Nt,1)=sparse_inverse;
%   for i=1:Nt
%     diag_block_inv_S(i).init((1/ctrl_inner11.omega)*diag_block_S{i},ctrl_inner11);
%   end

%   % set initial solution
%   x0=zeros(Np,1);
%   info_S=info_solver;

%   approach=controls.extra_info;
% 			       %approach="stationary_block_triangular"
%   if ( strcmp(approach, 'block_diag')) 
%     prec = @(x) apply_block_triangular_inverse(diag_block_inv_S,below_diag_block_S,'L',x);
%   elseif ( strcmp(approach, 'block_triangular') )
%     prec = @(x) apply_block_diag_inverse(diag_block_inv_S,x);
%   end

%   % solve 
%   tic;
%   [dp,info_S] = apply_iterative_solver(@(y) S*y, fp, ctrl_outer, prec,x0);
%   outer_cpu = toc;
  
  
%   % get solution of system JF d  = F  
%   dr = inv_C*(B2*dp-rhs(Np+1:Np+Nr));
%   ds = JF.ss\(h-JF.sr*dr);
%   d = [dp; dr; ds];

%   % store info
%   inner_nequ=diag_block_inv_S(1).nequ;
%   inner_iter=0;
%   for i  = 1:Nt
%     inner_iter=inner_iter+diag_block_inv_S(i).cumulative_iter;
%   end
%   outer_iter=info_S.iter

    
%   normd = norm(d);

elseif sol==10
  % Krylov method for full system
  % + 
  % use P=( SAC B1T ) as preconditoner 
  %       ( 0    -C )
  % cpu time preprocessing main solver
  time_prep=tic;

	% Ground matrix A, B1T and f1.
	% Change the upper bound cycle to see the effect of
	% "grounding" all blocks. 
	if (controls.ground == 1)
		for iblock = 1:1
			irow = controls.ground_node+(iblock-1)*ncellphi;
			A(irow,:) = sparse(Np,1);
			A(irow,irow)= 1.0;
			B1T_time(irow,:) = sparse(Nr,1);
			B1T_space(irow,:) = sparse(Nr,1);
			B1T(irow,:) = sparse(Nr,1);
			f1(irow)=0.0;
		end
	end

	W_mat=zeros(N,Nr);
	for k = 1:N
	  W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = JF.area2h';
  end
	

  % scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2}
  if (controls.diagonal_scaling)
		[sqrt_diag_A,sqrt_diag_C,inv_sqrt_diag_A,inv_sqrt_diag_C] = get_diagonal_scaling (A,C)
		A  =  diagA_scaling*A*  diagA_scaling;
		B1T= diagA_scaling*B1T*diagC_scaling;
		B1T_time  = diagA_scaling*B1T_time*diagC_scaling;
		B1T_space = diagA_scaling*B1T_space*diagC_scaling;
		B2 = diagC_scaling*B2* diagA_scaling;
		C  = diagC_scaling*C*  diagC_scaling;
		f1 = diagA_scaling * f1;
		f2 = diagC_scaling * f2;		
  end
  
	
  % init diag(C)^{-1}
	diag_C = spdiags(C,0);
  %fprintf('(%9.4e <= C <= %18.12e) \n',min(diag_C),max(diag_C))
  inv_diag_C = sparse(1:Nr,1:Nr,(1.0./diag_C)',Nr,Nr);
	% if C is diagonal we can invert easily. Otherwise we use agmg
	if ( nnz(C)==Nr)
		inverseC = @(x) inv_diag_C*x;
	else
		ctrl_loc=controls.ctrl_inner22;
		index_agmg=index_agmg+1;
		ctrl_loc.index_agmg=index_agmg;
		inverse_sparse_C=sparse_inverse;
    inverse_sparse_C.init(C+controls.relax22*speye(Nr,Nr),ctrl_loc);
		inverseC = @(x) inverse_sparse_C.apply(x);
	end

	P_mat=W_mat*inv_diag_C*B2;
	alphas=W_mat*inv_diag_C*f2;

		
	
	% Extra term B1T*inv_C*B2 of the primal schur complemenent
	% Form time-time, time-space and space-space component
	Stt = B1T_time* (inv_diag_C*B1T_time' );
	Stx = B1T_time* (inv_diag_C*B1T_space');
	Sxx = B1T_space*(inv_diag_C*B1T_space');

	
	% select
	if ( strcmp(controls.assembly_S,'full') )
		%disp('here')
		S = Stt+Stx+Stx'+Sxx;
  elseif ( strcmp(controls.assembly_S,'tt') ) 
    S = Stt;
  elseif ( strcmp(controls.assembly_S,'tt_xx') ) 
    S = Stt+Sxx;
	elseif ( strcmp(controls.assembly_S,'tx_xt') )
    S = Stx+Stx';
  elseif ( strcmp(controls.assembly_S,'tt_tx_xt') )
    S = Stt+Stx+Stx';
	elseif ( strcmp(controls.assembly_S,'harmonic') )
		% assume that phi is harmonic and use chain rule on
		% div(g\grad phi)=g*div(\gradphi) + grad phi * grad rho=grad phi * grad rho
		% BT~B
		% div invC B~= B
		%rhos = Rst*RHt*It*rho_all;
		% Rst = repmat({JF.Rs},1,N+1);
    % Rst = blkdiag(Rst{:});
		rec_time = assembleRHt(N-1,ncellrho);
		nei=size(JF.grad,2);
		gradt = assembleGradt(N,ncellphi,nei,JF.div');

		
		% crhos=spdiags(inv_diag_C,0);
		% assume slack is zeros for t=0 t=1,
		%rhos_all=[F.rho_in;spdiags(Dr,0);F.rho_f];
		%s_all  = [zeros(ncellrho,1);spdiags(Ds,0);zeros(ncellrho,1)];
		%m_all  = [area2h; spdiags(Mrho);area2h];

		%c_all  = s_all.*m_all./rhos_all;
		
		% avg in time and project to edges h
		te = size(gradt,1);
		%inv_C_cell = spdiags(RHt*JF.It*c_all,0,Np,Np);

		% G from edge to cellphi
	
		%JF.gradphi;
		%G= Rst'*spdiags(JF.gradphi,0,te,te);

		%		S=-JF.B1T_spaceG*gradt*inv_C*G*gradt;
		Rst = repmat({JF.Rs},1,N+1);
    Rst = blkdiag(Rst{:});
		I_cellphi_cellrho = assembleIt(N-1,ncellphi,ncellrho,JF.I);
		%  (N+1)*Kh<-(N+1)*e<-(N+1)*e<-(N+1)*Kh <-(N+1)*K2h <- N*K2h
		new_B1T_space= Rst' * spdiags(JF.gradphi,0,te,te) * gradt * I_cellphi_cellrho * rec_time';

		if (0)
			nei=size(JF.grad,2);
			gradt = assembleGradt(N-1,ncellphi,nei,JF.grad);
			Nr
			size(gradt)
			size(JF.G')
			size(rec_time')
			size(I_cellphi_cellrho)
			rec_time = assembleRHt(N-2,size(JF.grad,1));
			new_B1T_space= - I_cellphi_cellrho *  JF.G'*rec_time'* gradt ;
		end

		Stx = B1T_time* (inv_diag_C*B1T_space');
		Sxt = new_B1T_space*(inv_diag_C*B1T_time');
		Sxx = new_B1T_space*(inv_diag_C*B1T_space');

		S=Stt+Stx+Sxt+Sxx;

	elseif ( strcmp(controls.assembly_S,'tt_xx_lamped') ) 
    S = Stt
		% using a sort of mass lamping since weight in the mxx matrix
		% have the following structure
		% g1**2 g1*g2
		% g1*g2 2*g2**2 g1*g3
		%       
		i=1;
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)=...
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i  )*ncellphi+1:(i+1)*ncellphi,(i  )*ncellphi+1:(i+1)*ncellphi);
		% 
		for i=2:Nt-1
			S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)=...
			S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
			Sxx((i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
			Sxx((i-2)*ncellphi+1:(i-1)*ncellphi,(i-2)*ncellphi+1:(i-1)*ncellphi)+...
			Sxx((i  )*ncellphi+1:(i+1)*ncellphi,(i  )*ncellphi+1:(i+1)*ncellphi);
		end
		%
		i=Nt;
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)=...
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i-2)*ncellphi+1:(i-1)*ncellphi,(i-2)*ncellphi+1:(i-1)*ncellphi);

	elseif ( strcmp(controls.assembly_S,'blockdiagMxx') )
    block=sparse(Np,Np);
    for i=1:Nt
			block(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi)=...
			Sxx(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi);
    end
    S = block;
		
  elseif ( strcmp(controls.assembly_S,'tt_tx_xt_blockdiag_xx') )
    block=sparse(Np,Np);
    for i=1:Nt
			block(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi)=...
			Sxx(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi);
    end
    S = Stt+Stx+Stx'+block;
  end 

  % form Primal Schur complement (possibly approximated)
	S=A+S;


	% permute indeces primal block tominimize bandwidth
  if ( controls.permute == 1)
    % define permutation
    perm = symrcm(S);
    iperm(perm)=(1:Np)';
    
    A=A(perm,perm);
    B1T=B1T(perm,:);
		B1T_time=B1T_time(perm,:);
		B1T_space=B1T_space(perm,:);
    B2=B2(:,perm);
		f1=f1(perm);
		S=S(perm,perm);
  end
  

  % define (approximate) inverse of SAC=A+B1T diag(C)^{-1} B2 
  approach=controls.inverse11;  
  relax4prec=controls.relax_inv11;
	ctrl_inner11=controls.ctrl_inner11;

  timeS=tic;
  if (strcmp(approach,'full'))
    % init inverse of full S
		%disp('here')
		Nsys=1;
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
    inv_S.dimblock=ncellphi;
		%inv_S.info();
    inverse_cpu=inv_S.init_cpu;
		
    inverse_block11 = @(x) inv_S.apply(x);
    inner_nequ=Np;
  elseif(strcmp(approach,'block_triang'))
    %
    % partion diagonal and upper diagonal blocks
    %
		Nsys=Nt;
    diag_block_S       = cell(Nt, 1);
    below_diag_block_S = cell(Nt-1, 1);
    
    nAi=ncellphi;
    for i=1:Nt-1
      diag_block_S{i}      =S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
      below_diag_block_S{i}=S(i    *nAi+1 : (i+1)*nAi , (i-1)*nAi+1 : i*nAi);
    end
    i=Nt;
    diag_block_S{i}=S((i-1)*nAi+1:i*nAi,(i-1)*nAi+1:i*nAi);

    %				
    % define inverse on diagonal block
    %
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
    
    inverse_block11 = @(x) apply_block_triangular_inverse(diag_block_invS,below_diag_block_S,'L',x);
    inner_nequ=nAi;
  elseif(strcmp(approach,'block_diag'))
    nAi=ncellphi;
    inner_nequ=nAi;
    Nsys=Nt;
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
      diag_block_inv_S(i).init(S((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)+...
			      relax4prec*speye(nAi,nAi),ctrl_loc);
      diag_block_inv_S(i).name=sprintf('inverse_SAC%d%d',ctrl_inner11.label,i,i);
    end

    
    inverse_block11  = @(x) apply_block_diag_inverse(diag_block_inv_S,x);
  end
	preprocess_S=toc(timeS);

	prec_cpu=preprocess_S;
  preprocess_cpu=toc(time_prep);
	
  if( controls.solve_primal)
		if ( strcmp(controls.assembly_S,'full') & ...
				 controls.ctrl_inner11.tolerance<= controls.ctrl_outer.tolerance)
			
			fp = f1 + B1T*(inverseC(f2));
			sol_phi = inverse_block11(fp);
			sol_rho = inverseC(B2*sol_phi-f2);
			d=[sol_phi;sol_rho];
			if (controls.permute==1)
				d(1:Np) = d(iperm);
			end
			
			% get the ds 
			ds = Dr\(h-Ds*d(Np+1:Np+Nr));
			d=[d;ds];
		else
			disp('It not possible to directly solve the fully reduced system')
		end
  else
		%
		% Define action of preconditoner
		% 
		prec = @(v) SchurAC_based_preconditioner(v, inverse_block11, inverseC ,...
																						 @(y) B1T*y,@(x) B2*x,Np,Nr,controls.outer_prec,ncellphi);


		% solve with fgmres or bicgstab if prec is linear
		outer_timing=tic;
		rhs=[f1;f2];
		[d,info_J]=apply_iterative_solver(@(x) apply_saddle_point(x,@(y) A*y ,...
																															@(y) B1T*y,...
																															@(z) B2*z,...
																															@(z) C*z,...
																															Np,Nr),...
																			rhs, controls.ctrl_outer, prec );
		outer_iter=uint64(info_J.iter);
		outer_cpu=toc(outer_timing);		
	
		% permute back 
		if (controls.permute==1)
			d(1:Np) = d(iperm);
		end

		% rescale solutions
		if (controls.diagonal_scaling)
			d(1:Np)      = inv_sqrt_diag_A* d(1:Np);
			d(1+Np:Np+Nr)= inv_sqrt_diag_C * d(1+Np:Np+Nr);
		end


		% get the ds 
		ds = Dr\(h-Ds*d(Np+1:Np+Nr));
		d=[d;ds];

		
		
		% get inner solver info and free memory
		if(strcmp(approach,'full'))
			inner_iter=inv_S.cumulative_iter;
			inv_S.kill();
		elseif(strcmp(approach,'block_diag'))
			inner_iter=0;
			for i=1:Nt
				inner_iter = inner_iter + diag_block_invS(i).cumulative_iter;
				diag_block_inv_S(i).kill()
			end
		elseif(strcmp(approach,'block_triang'))
			inner_iter=0;
			for i=1:Nt
				inner_iter = inner_iter + diag_block_invS(i).cumulative_iter;
				diag_block_inv_S(i).kill()
			end
		end
		

		info_J.flag=1;
		if (info_J.res <= controls.ctrl_outer.tolerance)
			info_J.flag=0;
		end


		% get info 
		flag=info_J.flag;
		relres=info_J.res;
		iter=info_J.iter;	
		normd = norm(d);
		[ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);

		%print resume message
		resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.2e | nequ=%d - avg inner=%d cpu build S=%1.2e',...
											 info_J.flag,relres,resp,resr,ress,outer_iter,outer_cpu,...
											 inner_nequ, inner_iter/(outer_iter*Nsys), ...
											 preprocess_S);

	end

	
elseif sol==11
	% SIMPLE preconditioenr approach
	% Krylov method for full system
  % + 
  % use P=( I -(~A)^{-1} B2) ( (~A)^{-1}          ) ( I              )
  %       ( 0 I            ) (           (~S)^{-1}) (-(~A)^{-1} B1T I)
  %
  % 
  % based on
  % ~A=diag(A)
  % ~A=approximate inverse of A given by inner_solver
  % S=-C-B2*(diag(A))^{-1} B1T 


	% Ground matrix A, B1T and f1.
	% Change the upper bound cycle to see the effect of
	% "grounding" all blocks. 
	if (controls.ground == 1)
		for iblock = 1:1
			irow = controls.ground_node+(iblock-1)*ncellphi;
			A(irow,:) = sparse(Np,0);
			A(irow,irow)= 1.0;
			B1T_time(irow,:) = sparse(Nr,0);
			B1T_space(irow,:) = sparse(Nr,0);
			B1T(irow,:) = sparse(Nr,0);
			f1(irow)=0.0
		end
	end

  
	% scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
		[sqrt_diag_A,sqrt_diag_C,inv_sqrt_diag_A,inv_sqrt_diag_C] = get_diagonal_scaling (A,C)
		A  =  diagA_scaling*A*  diagA_scaling;
		B1T= diagA_scaling*B1T*diagC_scaling;
		B1T_time  = diagA_scaling*B1T_time*diagC_scaling;
		B1T_space = diagA_scaling*B1T_space*diagC_scaling;
		B2 = diagC_scaling*B2* diagA_scaling;
		C  = diagC_scaling*C*  diagC_scaling;
		f1 = diagA_scaling * f1;
		f2 = diagC_scaling * f2;		
  end
  
	% measure cpu time for prec assembly
  time_prec=tic;
	
  % define (~A)^{-1}  as diag(A)^{-1}
	inv_diag_A = 1.0./(spdiags(A,0)+controls.relax_inv11);
  inv_Diag_A = sparse(1:Np,1:Np,(inv_diag_A)',Np,Np);
	inv_A=@(x) inv_Diag_A*x;


	if (controls.study_eigen)

		for i=1:N
			nAi=ncellphi;
			matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi);
			inv_Diag_Ai = sparse(1:nAi,1:nAi,...
													 (inv_diag_A((i-1)*nAi+1 :     i*nAi))',...
													 nAi,nAi);
			eigenvalues=study_eigenvalues(matrixAi*inv_Diag_Ai, '(S)^{-1}droped',1);
			fprintf('eig(2)=%1.2e, eig(max)=%1.2e cond=%1.2e\n',...
							eigenvalues(2),eigenvalues(nAi),eigenvalues(nAi)/eigenvalues(2))
		end
	end

	
 
  % assembly (~S)^{-1} = ( -(C+B2 * diag(A)^{-1} B1T) )^{-1} 
	build_S=tic;

  SCA = (C+B2*inv_Diag_A*B1T);

  inv_SCA=sparse_inverse;
  inv_SCA.init(SCA+controls.relax_inv22*speye(Nr,Nr),controls.ctrl_inner22);
  inv_SCA.info_inverse.label='inv_tildeS';
  inv_SCA.cumulative_iter=0;
  inv_SCA.cumulative_cpu=0;
  inverse_block22 = @(x) -inv_SCA.apply(x);
  cpu_assembly_inverseS=toc(build_S);
	  
	
  % store time required to build prec
  prec_cpu=toc(time_prec);
	
  % set dimension of inner solver
  inner_nequ=Nr;

	
	% Define action of block preconditoner 
	prec = @(x) SchurCA_based_preconditioner(x, inv_A,...
																					 inverse_block22,...
																					 @(y) B1T*y,...
																					 @(z) B2*z,...
																					 Np,Nr,...
																					 controls.outer_prec,ncellphi);

	% call iterative solver
	outer_timing=tic;
	rhs=[f1;f2];

	[d,info_J]=apply_iterative_solver(@(x) apply_saddle_point(x,@(y) A*y ,...
																														@(y) B1T*y,...
																														@(z) B2*z,...
																														@(z) C*z,...
																														Np,Nr), ...
																		rhs, controls.ctrl_outer, @(z) prec(z),[],controls.left_right );
	outer_cpu=toc(outer_timing);
	
	if (controls.diagonal_scaling)
		d(1:Np)      = inv_diag_A*d(1:Np);
		d(1+Np:Np+Nr)= inv_diag_C*d(1+Np:Np+Nr);  
	end

	% get the ds 
	ds = Dr\(h-Ds*d(Np+1:Np+Nr));
	d=[d;ds];
  

	% get info
	inner_iter=inv_SCA.cumulative_iter;
	outer_iter=uint32(info_J.iter);
	
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  relres=res;
  normd = norm(d);
  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| S : nequ=%d inn=%d ass=%1.2e',...
   		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
   		     inv_SCA.nequ, inv_SCA.cumulative_iter,cpu_assembly_inverseS);

	% free memory of inner solvers
  inv_SCA.kill();
elseif sol==110
	% SIMPLE preconditioenr to the full system 
	
  % copy controls
  indc=controls.indc;
  sol=controls.sol;
  ctrl_inner22=controls.ctrl_inner22;
  ctrl_outer=controls.ctrl_outer;
  verbose=controls.verbose;
 
 
	% scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  %if (controls.diagonal_scaling)
  %  [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  %end
  

  time_prec=tic;
	
	% define approxiamte inverse of (~A)^{-1}  
	inv_diagA = 1.0./(spdiags(A,0)+controls.relax4inv11);
  inv_DiagA = sparse(1:Np,1:Np,(invdiagA)',Np,Np);
	inv_A=@(x) invDiagA*x;

 
  
  % assembly S22=(B2 * diag(A)^{-1} B1T) 
  S22 = (B2*inv_DiagA*B1T);
	build_S=tic;
  % init inverse of full S
  if ( verbose >= 2)
    fprintf('INIT inverse S\n')
  end
  inv_S22=sparse_inverse;
  inv_S22.init(S22+controls.relax4inv22*speye(Nr,Nr),ctrl_inner22);
  inv_S22.info_inverse.label='schur_ca';
  inv_S22.cumulative_iter=0;
  inv_S22.cumulative_cpu=0;
  inverse_block22 = @(x) -inv_S22.apply(x);

	
	% S33^{-1}=(Dr                  +Ds*S^{-1}*M)^{-1}
	%         =(Dr*M^{-1}*S*S^{-1}*M+Ds*S^{-1}*M)^{-1}
	%         =M^{-1} S * (Dr*M^{-1}*S+Ds)^{-1}
	inv_M = sparse(1:Nr,1:Nr,(1./spdiags(M,0))',Nr,Nr);			
	S33  = (Dr*inv_M*S22 + Ds);

	
	inv_S33=sparse_inverse;
  inv_S33.init(S33+controls.relax4inv22*speye(Nr,Nr),ctrl_inner22);
  inv_S33.info_inverse.label='schur_ca';
  inv_S3.cumulative_iter=0;
  inv_S33.cumulative_cpu=0;
  if ( verbose >= 2)
    fprintf('DONE inverse S\n') 
  end 
  inverse_block33 = @(x) inv_M*(S22*(inv_S33.apply(x)));
  cpu_assembly_inverseS=toc(build_S);
	
  % store time required to build prec
  prec_cpu=toc(time_prec);
	
  % set dimension of inner solver
  inner_nequ=Nr;
	inner_nequ2=Nr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % END DUAL SCHUR COMPLEMENT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Define action of block preconditoner 
	prec11_22 = @(x) SchurCA_based_preconditioner(x, inv_A,...
																								inverse_block22,...
																								@(y) B1T*y,...
																								@(z) B2*z,...
																								Np,Nr,...
																								'full',ncellphi);

	
	% global operator
	global_A = @(x) apply_saddle_point(x,@(y) A*y ,...
																				 @(y) B1T*(y),...
																				 @(z) (B2*z),...
																				 @(z) -(R*(z)),...
																				 Np,Nr);
	global_B1T = @(y) [zeros(Np,1);(M*y)];
	global_B2  = @(x) Ds*(x(Np+1:Np+Nr));
	global_C   = @(y) -Dr*y;

	J=@(x) apply_saddle_point(x,global_A ,...
														global_B1T,...
														global_B2,...
														global_C,...
														Np+Nr,Nr);
	
	% Define action of block preconditoner 
	prec = @(x) SchurCA_based_preconditioner(x, prec11_22 ,...
																					 @(y) inverse_block33(y),...
																					 global_B1T,...
																					 global_B2,...
																					 Np+Nr,Nr,...
																					 controls.outer_prec,ncellphi);
	
	
	outer_timing=tic;
	rhs=[f;g;h];
	[d,info_J]=apply_iterative_solver(J, rhs, ctrl_outer, @(z) prec(z),[],controls.left_right );
	outer_cpu=toc(outer_timing);
	
	
	% get info
	inner_iter=inv_S22.cumulative_iter;
	outer_iter=uint32(info_J.iter);

	inner_iter2=inv_S33.cumulative_iter;
  
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;


 

  relres=res;
  normd = norm(d);
  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| S : nequ=%d inn=%d ass=%1.2e',...
   		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
   		     Nr, inner_nequ,cpu_assembly_inverseS);

  inv_S22.kill();
	inv_S33.kill();

	
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
    diag_block_inv_A(Nt,1)=sparse_inverse;

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
      diag_block_inv_A(i).name=sprintf('inverse A%d',i);
      
      % get block add relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi)+relax4prec*speye(nAi,nAi);
      
				% define inverse
      diag_block_inv_A(i).init( matrixAi,ctrl_loc);     

    end
% define function A shifter( inverse A grounded ( rhs set to zero in some nodes )
    inverse11  = @(x) apply_inverseAtilde(diag_block_inv_A,A,vectors_x,indeces_global,indeces_local,x);
  end

  cpu_assembly_inverseA=toc(build_inverseA);
  preprocess_cpu=toc;

  %
				%
  build_inverseS=tic;
  approach_schurCA=controls.approach_schurCA;
  if (strcmp(approach_schurCA,'diagA'))
    % define explicitely inverse of diag(A)
    inv_DiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);

    % assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
    approx_SCA = (C+B2*inv_DiagA*B1T);

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
    inv_DiagA = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);

    % assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
    approx_SCA = (C+B2*inv_DiagA*B1T);

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
    inv_DiagA = sparse(1:Np,1:Np,(1.0./DiagA)',Np,Np);
	
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
    
    if (swap_sign)
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
	% The system is partially reduced to a standard saddle point system.
	% The approach from : "A PRECONDITIONER FOR GENERALIZED SADDLE POINT PROBLEMS"
  % solve M=(A B1T) = H + S
  %         (-B2 C) 
  % where 
  % H=(A  ) S=(     B1T)
  %   (  C)   (-B2     )
  % Krylov method for full system
  % + 
  % use P = 1/(2a)*(H+aI)(S+aI) or some variations
  %
	
  % cpu time preprocessing main solver
  tic;

	% handle singularity of the system
	if (controls.ground == 1)
		% Ground matrix A, B1T and f1.
		% Change the upper bound cycle to see the effect of
		% "grounding" all blocks. 
		for iblock = 1:1
			irow = controls.ground_node+(i-1)*ncellphi;
			A(irow,:) = sparse(Np,0);
			A(irow,irow)= 1.0;
			B1T(irow,:) = sparse(Nr,0);
			f(irow) = 0;
		end
	end

	
  time_manipulate=tic;
  
  % scale system by diag(A,C)^{-1/2} M  diag(M)^{-1/2} diag(A,C)^{1/2}
	% x = diag(A,C)^{-1/2} rhs
  if (controls.diagonal_scaling)
		[sqrt_diag_A,sqrt_diag_C,inv_sqrt_diag_A,inv_sqrt_diag_C] = get_diagonal_scaling (A,C);
		A  = inv_sqrt_diag_A*A*  inv_sqrt_diag_A;
		B1T= inv_sqrt_diag_A*B1T*inv_sqrt_diag_C;
		B1T_time  = inv_sqrt_diag_A*B1T_time*inv_sqrt_diag_C;
		B1T_space = inv_sqrt_diag_A*B1T_space*inv_sqrt_diag_C;
		B2 = inv_sqrt_diag_C* B2* inv_sqrt_diag_A;
		C  = inv_sqrt_diag_C *C * inv_sqrt_diag_C;
		f1 = inv_sqrt_diag_A * f1;
		f2 = inv_sqrt_diag_C * f2;
  end
  cpu_manipulate=toc(time_manipulate);
  
  % change sign of rhs
  %matrixM = [A B1T; -B2 C];
	
  %get relazation paramter
  alpha=controls.alpha;

  time_prec=tic;
  timeA=tic;
  % define inverse of H+alpha*I

	% define inverse of A
	% the 'full' approach is left for debugging purpose
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
	beta=alpha;
  inv_Calpha=sparse_inverse;
  inv_Calpha.init( C + beta*speye(Nr,Nr) , controls.ctrl_innerC);

  % concatenate
  inv_Halpha = @(x) [inv_Aalpha(x(1:Np)); inv_Calpha.apply(x(Np+1:Np+Nr))];

  % Define inverse of S_alpha=S+(alpha*I;beta*I)=(alpha*I  B1T     )
	%                                              (-B2      -alpha*I)
  timeS=tic;
  if (strcmp(controls.approach_inverse_S,'primal'))
		beta
		if (beta > 1.0e-10)
			% invert S using its primal schur complement
			S_primal=alpha*speye(Np,Np)+1/beta*B1T*B2;
			inv_S11=sparse_inverse;
			inv_S11.init(S_primal,controls.ctrl_innerS);
			
			inv_Salpha = @(x) SchurAC_based_preconditioner(x,...
																										 @(xin) inv_S11.apply(xin),...
																										 @(y) -y/beta ,...
																										 @(y) B1T*y,...
																										 @(x)-B2*x,...
																										 Np,Nr,...
																										 'full',ncellphi);
		else
			disp('beta too small')
			return
		end
	elseif (strcmp(controls.approach_inverse_S,'dual'))
		% invert S using its dual schur complement
    S_dual=beta*speye(Nr,Nr)+1/alpha*B2*B1T;
    inv_S11=sparse_inverse;
    inv_S11.init(S_dual,controls.ctrl_innerS);
    inv_Salpha = @(x) SchurCA_based_preconditioner(x, @(xin) xin/alpha,...
						   @(yin) inv_S11.apply(yin),...
						   @(y) B1T*y,...
						   @(z) -B2*z,...
						   Np,Nr,...
						   'full',ncellphi);
  end
  cpu_assembly_inverseS=toc(timeS);
  

	% define preconditioner if P=HS or P=SH
  if (strcmp(controls.approach_prec,'HS'))
    prec = @(x) inv_Salpha(inv_Halpha(x))/(2*alpha);
  elseif ( strcmp(controls.approach_prec,'SH'))
    prec = @(x) inv_Halpha(inv_Salpha(x))/(2*alpha);
  end
 
  %preprocess finished
  preprocess_cpu=toc;
  prec_cpu=toc(time_prec);

  % solve linear system with Krlyov method
  outer_timing=tic;
	rhs=[f1;-f2];
  [d,info_J]=apply_iterative_solver(@(v) apply_saddle_point(v,@(x) A*x ,...
																														@(y) B1T*y,...
																														@(x) -B2*x,...
																														@(y) -C*y,... % there is a minus inside
																														Np,Nr,N),...
																		rhs, controls.ctrl_outer, prec,[],controls.left_right );
	% scale back the system
  if (controls.diagonal_scaling)
    d(1:Np)      =inv_sqrt_diag_A*d(1:Np);
    d(1+Np:Np+Nr)=inv_sqrt_diag_C*d(1+Np:Np+Nr);  
  end
	% get the ds 
	ds = Dr\(h-Ds*d(Np+1:Np+Nr));
	d=[d;ds];
	
  outer_cpu=toc(outer_timing);

	% 
	% get statistics and free memory
	%
  if ( strcmp(approach_inverse_A,'full'))
    inner_iter_A=invAalpha.cumulative_iter;
    inner_nequ_A=invAalpha.nequ;
		invA.kill()
  elseif ( strcmp(approach_inverse_A,'block'))
    inner_nequ_A=diag_block_invA(1).nequ;
		inner_iter_A=0;
		inner_cpu_A=0;
		for i=1:Nt
      inner_iter_A = inner_iter_A + diag_block_invA(i).cumulative_iter;
			inner_cpu_A  = inner_cpu_A + diag_block_invA(i).cumulative_cpu;
			diag_block_invA(i).kill();
    end
  end
	inner_nequ=inner_nequ_A;
  inner_iter=inner_iter_A;

	
	inner_iter_S=inv_S11.cumulative_iter;
  inner_nequ_S=inv_S11.nequ;
  inner_nequ2=inner_nequ_S;
  inner_iter2=inner_iter_S;
	inv_S11.kill();	


  inner_nequ3=Nr;
  inner_iter3=inv_Calpha.cumulative_iter;
	inv_Calpha.kill();
	
  
  flag=info_J.flag;
  relres=info_J.res;
  outer_iter=info_J.iter;

  normd = norm(d);

  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);

	%
	% resume message
	%
  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| %s A: nequ=%d inn.=%d avg=%d ass=%1.2e | C : %s inn.=%d | S : nequ=%d inn=%d ass=%1.2e',...
		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
		     approach_inverse_A,inner_nequ_A,inner_iter_A,round(inner_iter_A/(Nt*outer_iter)),cpu_assembly_inverseA, ...
		     inv_Calpha.ctrl.label,inv_Calpha.cumulative_iter,...
		     inner_nequ_S, inner_iter_S,cpu_assembly_inverseS);
  
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

  % scale system by diag(M)^{-1/2} M  diag(M)^{-1/2} diag(M)^{1/2} x = diag(M)^{-1/2} rhs
  if (controls.diagonal_scaling)
    [A,B1T,B2,C,rhs,diagA_scaling,diagC_scaling]=scaling_system(A,B1T,B2,C,rhs);
  end

	
  if (strcmp(controls.gamma,'auto'))
    gamma=norm(A,'fro')/(norm(B1T,'fro')^2);
  else
    gamma=controls.gamma;
  end

	
	% set diagonal matrix W and its inverse
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
  end
    

  build_inv_augA=tic;
	% Extra term B1T*invM*B2 of the primal schur complemenent
	% Form time-time, time-space and space-space component
	Stt = JF.B1T_time* (invW*JF.B1T_time' );
	Stx = JF.B1T_time* (invW*JF.B1T_space');
	Sxx = JF.B1T_space*(invW*JF.B1T_time');

	% sel
	if ( strcmp(controls.assembly_S,'full') )
		S = Stt+Stx+Stx'+Sxx
  elseif ( strcmp(controls.assembly_S,'tt') ) 
    S = Stt;
  elseif ( strcmp(controls.assembly_S,'tt_xx') ) 
    S = Stt+Sxx;
	elseif ( strcmp(controls.assembly_S,'tx_xt') )
    S = Stx+Stx';
  elseif ( strcmp(controls.assembly_S,'tt_tx_xt') )
    S = Stt+Stx+Stx';
	elseif ( strcmp(controls.assembly_S,'tt_xx_lamped') ) 
    S = Stt
		% using a sort of mass lamping since weight in the mxx matrix
		% have the following structure
		% g1**2 g1*g2
		% g1*g2 2*g2**2 g1*g3
		%       
		i=1;
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)=...
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i  )*ncellphi+1:(i+1)*ncellphi,(i  )*ncellphi+1:(i+1)*ncellphi);
		% 
		for i=2:Nt-1
			S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)=...
			S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
			Sxx((i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
			Sxx((i-2)*ncellphi+1:(i-1)*ncellphi,(i-2)*ncellphi+1:(i-1)*ncellphi)+...
			Sxx((i  )*ncellphi+1:(i+1)*ncellphi,(i  )*ncellphi+1:(i+1)*ncellphi);
		end
		%
		i=Nt;
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)=...
		S(    (i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i-1)*ncellphi+1:(i  )*ncellphi,(i-1)*ncellphi+1:(i  )*ncellphi)+...
		Sxx((i-2)*ncellphi+1:(i-1)*ncellphi,(i-2)*ncellphi+1:(i-1)*ncellphi);

	elseif ( strcmp(controls.assembly_S,'blockdiagMxx') )
    block=sparse(Np,Np);
    for i=1:Nt
			block(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi)=...
			Sxx(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi);
    end
    S = block;
		
  elseif ( strcmp(controls.assembly_S,'tt_tx_xt_blockdiag_xx') )
    block=sparse(Np,Np);
    for i=1:Nt
			block(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi)=...
			Sxx(1+(i-1)*ncellphi:i*ncellphi,1+(i-1)*ncellphi:i*ncellphi);
    end
    S = Stt+Stx+Stx'+block;
  end 


	% form matrix A + gamma S
	% 
  augA=A+gamma*S;
 

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
		%inv_augA.ctrl.info()
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
  prec = @(x) SchurCA_based_preconditioner(x, ...
																					 inverse11,...
																					 inverse22,...
																					 @(y) augB1T*y,...
																					 @(z) B2*z,...
																					 Np,Nr,...
																					 controls.outer_prec,ncellphi);



  preprocess_cpu=toc(preprocess);

  
  % solve with fgmres 
  outer_timing=tic;
  

  
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

    
		%augjacobian = [augA augB1T; B2 -C];
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
  [d,info_J]=apply_iterative_solver(@(v) ...
																		 apply_saddle_point(v,@(x) augA*x ,...
																												@(y) augB1T*y,...
																												@(x) B2*x,...
																												@(y) C*y,...
																												Np,Nr,N), ...
																		augrhs, ctrl_outer, prec );
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
	W_mat=zeros(N,Nr);
	for k = 1:N
    W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = area2h';
  end
	PT = @(y) PT_apply(y,area2h);

	
	diag_C = spdiags(C,0);
	inv_diag_C = sparse(1:Nr,1:Nr,(1.0./diag_C)',Nr,Nr);
	P_mat=W_mat*inv_diag_C*B2;
	alphas=W_mat*inv_diag_C*f2;

	P_perturbation=sparse(Np,Np);
	f_perturbation=zeros(Np,1);
	for i=1:N
		irow=(i-1)*ncellphi+1;
		P_perturbation(irow,:)= P_mat(i,:);
		f_perturbation(irow)  = alphas(i);
	end
	i=N+1;
	irow=(i-1)*ncellphi+1;
	P_perturbation(irow,Np-ncellphi+1:Np)=P_mat(N,Np-ncellphi+1:Np);
	f_perturbation(irow)  = 0;
	
	inv_diag_A = 1.0./(spdiags(A,0)+controls.relax_inv11);
  inv_Diag_A = sparse(1:Np,1:Np,(inv_diag_A)',Np,Np);
	% add perturbation
	A=A+P_perturbation;
	f1=f1+f_perturbation;


	time_prec=tic;
	
	% partion matrix 
  nAi=ncellphi;
  
  % use block inverse
  diag_block_invA(Nt,1)=sparse_inverse;
  ctrl_loc=ctrl_solver; 
  for i=1:Nt
		ctrl_inner11=controls.ctrl_inner11;
		index_agmg=index_agmg+1;
    ctrl_loc=ctrl_solver;
    ctrl_loc.init(ctrl_inner11.approach,...
    							ctrl_inner11.tolerance,...
    							ctrl_inner11.itermax,...
    							ctrl_inner11.omega,...
    							ctrl_inner11.verbose,...
    							sprintf('%sA%d',ctrl_inner11.label,i),...
									index_agmg);
    diag_block_invA(i).name=sprintf('inverse A%d',i);

    % create local block and passing to solver
    % with a potential relaxation
    matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + ...
						 controls.relax_inv11*speye(nAi,nAi) ;
    diag_block_invA(i).init(matrixAi,ctrl_loc);

		nsysA=Nt;
  end

	extra_diagonal=cell(N,1);
	for i=1:N
		extra_diagonal_matrices{i}=A((i-1)*nAi+1 :     i*nAi , (i)*nAi+1 : (i+1)*nAi);
	end
		
  % define function
	inv_A  = @(x) apply_block_triangular_inverse(diag_block_invA,extra_diagonal_matrices,'U',x);
	
	% assembly (~S)^{-1} = ( -(C+B2 * diag(A)^{-1} B1T) )^{-1} 
	build_S=tic;

	%inv_diag_A = 1.0./(spdiags(A,0)+controls.relax_inv11);
  %inv_Diag_A = sparse(1:Np,1:Np,(inv_diag_A)',Np,Np);
  SCA = (C+B2*inv_Diag_A*B1T);

  inv_SCA=sparse_inverse;
  inv_SCA.init(SCA+controls.relax_inv22*speye(Nr,Nr),controls.ctrl_inner22);
  inv_SCA.info_inverse.label='inv_tildeS';
  inv_SCA.cumulative_iter=0;
  inv_SCA.cumulative_cpu=0;
  inverse_block22 = @(x) -inv_SCA.apply(x);
  cpu_assembly_inverseS=toc(build_S);
	  
	
  % store time required to build prec
  prec_cpu=toc(time_prec);
	
  % set dimension of inner solver
  inner_nequ=Nr;

	
	% Define action of block preconditoner 
	prec = @(x) SchurCA_based_preconditioner(x, inv_A,...
																					 @(y) PT(inverse_block22(y)),...
																					 @(y) B1T*PT(y),...
																					 @(z) B2*z,...
																					 Np,Nr,...
																					 controls.outer_prec,ncellphi);

	% call iterative solver
	outer_timing=tic;
	rhs=[f1;f2];
	[d,info_J]=...
	apply_iterative_solver(@(x) apply_saddle_point(x,@(y) A*y ,...
																								 @(y) B1T*(y),...
																								 @(z) B2*z,...
																								 @(z) C*(z),...
																								 Np,Nr), ...
												 rhs, ...
												 controls.ctrl_outer, ...
												 @(z) prec(z),[],controls.left_right );
	outer_cpu=toc(outer_timing);
	

	% get the ds 
	ds = Dr\(h-Ds*d(Np+1:Np+Nr));
	d=[d;ds];
  

	% get info
	inner_iter=inv_SCA.cumulative_iter;
	outer_iter=uint32(info_J.iter);
	
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;

  relres=res;
  normd = norm(d);
  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| S : nequ=%d inn=%d ass=%1.2e',...
   		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
   		     inv_SCA.nequ, inv_SCA.cumulative_iter,cpu_assembly_inverseS);

	% free memory of inner solvers
  inv_SCA.kill();


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

    eigenvalues=study_eigenvalues(prec_matrix,'~S^{-1} S',1);

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
		f_lambda(k) = v' * get_slice(rho,k, ncellrho, N) - 1; % F^k_lambda
    W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = v'; % d F^k_lambda drho_k 
  end
  
  % build projector matrix
	factor=1/(norm(v))^2;
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
	Sys(indc,:) = sparse(1,dim); Sys(:,indc) = sparse(dim,1); Sys(indc,indc) = 1; rhs(indc)=0; 
	increment=Sys\rhs;
  
  % get phi,rho,s increment and lambda
  d=increment(1:Np+Nr+Nr);
  dl=increment(1+Np+Nr+Nr: Np+Nr+Nr+N);
	[resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
	resume_msg=sprintf('res=%1.1e | p=%1.1e r=%1.1e s=%1.1e\n',resnorm,resp,resr,ress);

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

	% assembly lambda residuum
	% f(k)=(area*rho_k-1)
	v = JF.area2h;
	f_lambda = zeros(N,1);
	W_mat=zeros(N,Nr);
	for k = 1:N
		f_lambda(k) = v' * get_slice(rho,k, ncellrho, N) - 1; % F^k_lambda
    W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = v'; % d F^k_lambda drho_k 
  end
  
  % build projector matrix
	factor=1/(norm(v))^2;
  P = @(y) P_apply(y,area2h);
	PT = @(y) PT_apply(y,area2h);

	ground_A=1;
	if ground_A==1
		%diagA=spdiags(SysR(1:Np,1:Np),0);
		node=1;
		dim=Np+Nr;
		for i=1
			%[dmax,imax]=max(diagA(1+(i-1)*ncellphi:i*ncellphi));
			%indc=imax+(i-1)*ncellphi;
			indc=node+(i-1)*ncellphi;
			
			A(indc,:) = sparse(1,Np);
			%A(:,indc) = sparse(dim,1);
			A(indc,indc) = 1;
			B1T(indc,:) = sparse(1,Nr);
			
			f1(indc)=0;		
		end
	end
	

	
	
	exact=1;
	if (exact ==1)
		SysR = sparse([...
  									A,   B1T;
  									B2, -C;   ...
  								]);
		% 
		relax=0e-12;
		SysR(1:Np,1:Np) = SysR(1:Np,1:Np) + relax*speye(Np);
		SysR(1+Np:Np+Nr,1+Np:Np+Nr) = SysR(1+Np:Np+Nr,1+Np:Np+Nr) - relax*speye(Nr);

		rhsR = [f1; f2];
		d_exact = SysR\rhsR;
		ds = Dr\(h-Ds*d_exact(Np+1:Np+Nr));
		d_exact=[d_exact;ds];
	end

	% defined H
	H=speye(ncellrho-1,ncellrho);   
	H(1:ncellrho-1,ncellrho)=-area2h(1:ncellrho-1)/area2h(ncellrho);
	mat_H = repmat({H},1,N);
	mat_H = blkdiag(mat_H{:});

	% defined k
	K=speye(ncellrho-1,ncellrho);   
	mat_K = repmat({K},1,N);
	mat_K = blkdiag(mat_K{:});
	%mat_K (ncellrho-1,ncellrho)=1.0;
	mat_K=mat_H;
	
	SysR = sparse([...
  										A,   B1T*mat_H';
  										mat_K*B2, -mat_K*C*mat_H';   ...
  									]);

	rhsR=blkdiag(speye(Np,Np),mat_K)*rhsR;
	d=SysR\rhsR;
	
	d = blkdiag(speye(Np,Np),mat_H')*d;
	ds = Dr\(h-Ds*d(Np+1:Np+Nr));
	d=[d;ds];

	figure
	plot(d_exact-d)


	% get lambda increment dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)
	dl = get_dl(JF,F,d);
  
  % correction of phi increment 
	d(1:Np) = dp_correction(d(1:Np),JF,dl);


	figure
	
	size(d_exact)
	size(d)
	plot(d_exact-d)

	figure
	spy(SysR)
	
	
	return
	

	

	mode=1;
	if (mode==0)
		rhsR = -[F.p; P*F.r; F.s];
		SysR = sparse([...
  									JF.pp,         JF.pr,   JF.ps, ;
  									P*JF.rp,       P*JF.rr, P*JF.rs;
  									sparse(Nr,Np), JF.sr,  JF.ss   ...
  								]);
		
		
		d = SysR\rhsR;
		res=norm(SysR*d - rhsR);
		fprintf('res aug=%1.2e\n',res)

	elseif(mode==1)
		%print_imbalance(f1,ncellphi)
		factors=zeros(N,1);
		for k=1:N
			factors(k)=1./norm(W_mat(k,:))^2;
		end
		P=speye(Nr,Nr)-W_mat'*diag(factors)*W_mat;
		
		% scale system
		diagonal_scaling=0;
		if diagonal_scaling==1
			[sqrt_diagA,sqrt_diagC,inv_sqrt_diagA,inv_sqrt_diagC] = get_diagonal_scaling (A,C);
			[A,B1T,B2,C,f1,f2,W_mat]=scale_saddle_point(A,B1T,B2,C,f1,f2,W_mat,inv_sqrt_diagA,inv_sqrt_diagC);
			for k=1:N
				factors(k)=1./norm(W_mat(k,:))^2;
			end
			P=speye(Nr,Nr)-W_mat'*diag(factors)*W_mat;
		end
		rhsR = [f1; P*f2];
		
		remove_singularity=3;
		if ( remove_singularity==1)
			disp('here')
			SysR = sparse([...
  										A,   B1T*P;
  										P*B2, -P*C*P;   ...
  									]);

			ground_A = 1;
			if ground_A==1
				%diagA=spdiags(SysR(1:Np,1:Np),0);
				node=1;
				dim=Np+Nr;
				for i=1:N+1
					%[dmax,imax]=max(diagA(1+(i-1)*ncellphi:i*ncellphi));
					%indc=imax+(i-1)*ncellphi;
					indc=node+(i-1)*ncellphi;
					
					SysR(indc,:) = sparse(1,dim);
					SysR(:,indc) = sparse(dim,1);
					SysR(indc,indc) = 1;
					
					rhsR(indc)=0;		
				end
			end
			
			% 
			relax=1e-10;
			%SysR(1:Np,1:Np) = SysR(1:Np,1:Np) + (1-ground_A)*relax*speye(Np);
			SysR(1+Np:Np+Nr,1+Np:Np+Nr) = SysR(1+Np:Np+Nr,1+Np:Np+Nr) - relax*speye(Nr);

			d = SysR\rhsR;
			
			
			
		elseif ( remove_singularity==2)
		
			SysR = sparse([...
  										A,   B1T;
  										B2, -C;   ...
  									]);

			change=[...
			 				 speye(Np), sparse(Np,Nr);...
			 				 sparse(Nr,Np), P; ....
			 			 ];

			
			node=1;
			dim=Np+Nr;
			for i=1:N+1
				indc=node+(i-1)*ncellphi;
				SysR(indc,:) = sparse(1,dim);
				SysR(:,indc) = sparse(dim,1);
				SysR(indc,indc) = 1;
				rhsR(indc)=0;
			end

			SysR = change*SysR;

			relax=1e-10;
			SysR(1+Np:Np+Nr,1+Np:Np+Nr) = SysR(1+Np:Np+Nr,1+Np:Np+Nr) - relax*speye(Nr);
						
			
			d = SysR\rhsR;

		elseif ( remove_singularity==3)
			H=speye(ncellrho-1,ncellrho);   
			H(1:ncellrho-1,ncellrho)=-area2h(1:ncellrho-1)/area2h(ncellrho);
			
			% H=block diagonal (P) 
			mat_H = repmat({H},1,N);
			mat_H = blkdiag(mat_H{:});

			SysR = sparse([...
  										A,   B1T*mat_H';
  										mat_H*B2, -mat_H*C*mat_H';   ...
  									]);

			rhsR=blkdiag(speye(Np,Np),mat_H)*rhsR;

			ground_A=0;
			if (ground_A)
				node=1;
				dim=size(SysR,1);
				for i=1:N+1
					indc=node+(i-1)*ncellphi;
					SysR(indc,:) = sparse(1,dim);
					SysR(:,indc) = sparse(dim,1);
					SysR(indc,indc) = 1;
					rhsR(indc)=0;
				end
			end
			relax=1e-10;
			SysR(1:Np,1:Np) = SysR(1:Np,1:Np) + (1-ground_A)*relax*speye(Np);
		

			
			d = SysR\rhsR;

			fprintf('res reduced=%1.1e \n',norm(SysR*d-rhsR));
			
			d = blkdiag(speye(Np,Np),mat_H')*d;
			norm(d(Np+1:Np+Nr)-ort_proj(d(Np+1:Np+Nr),W_mat))

		end
		
	
		if (diagonal_scaling == 1)
			d(1:Np)      = inv_sqrt_diagA * d(1:Np);
			d(1+Np:Np+Nr)= inv_sqrt_diagC * d(1+Np:Np+Nr);
			W_mat=W_mat*sqrt_diagC;
		end

		%figure
		%plot(d-d_exact(1:Np+Nr))
		

		ds = Dr\(h-Ds*d(Np+1:Np+Nr));
		d = [d; ds];
	elseif (mode==2)
		
		S=A+ B1T * ((P*C) \ (P*B2));
		fp=f1-B1T*((P*C)\f2);
		dp = S\fp;
		res=norm(S*dp - fp);
		fprintf('res aug=%1.2e\n',res)
		dr = (P*C)\(-f2-P*B2*dp);
		ds = Dr\(h-Ds*dr);
		d=[dp;dr;ds];
	end

	norm(d(Np+1:Np+Nr)-ort_proj(d(Np+1:Np+Nr),W_mat))
  d(Np+1:Np+Nr)=PT(d(Np+1:Np+Nr));
		
	% get lambda increment dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)
	dl=get_dl(JF,F,d);
  
  % correction of phi increment 
	d(1:Np) = dp_correction(d(1:Np),JF,dl);

	
	% compute new residuum
  [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);

	for i = 1:N
		imb=d(Np+1+(i-1)*ncellrho:Np+i*ncellrho)'*area2h;
		state_message=sprintf('y*area=%+1.1e',imb);
		%fprintf('%s \n',state_message);
	end

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

elseif (sol==19)
 	% get variable for assembling system
  rho=-spdiags(JF.ss);
  deltat=1/(N+1);

	% assembly lambda residuum
	% f(k)=(area*rho_k-1)
	v = JF.area2h;
	f_lambda = zeros(N,1);
	W_mat=zeros(N,Nr);
	for k = 1:N
		f_lambda(k) = v' * get_slice(rho,k, ncellrho, N) - 1; % F^k_lambda
    W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = v'; % d F^k_lambda drho_k 
  end
  
  rhsR = [f1; ort_proj(f2,W_mat)];

	ground_A = 1;
	if ground_A==1
		diagA=spdiags(A,0);
		node=1;
		dim=Np+Nr;
		for i=1:N+1
			[dmax,imax]=max(diagA(1+(i-1)*ncellphi:i*ncellphi));
			%indc=imax+(i-1)*ncellphi;
			indc=node+(i-1)*ncellphi;
			
			A(indc,:) = sparse(1,Np);
			A(:,indc) = sparse(Np,1);
			A(indc,indc) = 1;

			B1T(indc,:) = sparse(1,Nr);
			B2(:,indc) = sparse(Nr,1);
			
			rhsR(indc)=0;		
		end
	end

	diagonal_scaling=1;
	if diagonal_scaling==1
		inv_sqrt_A   = sparse(1:Np,1:Np,(1.0./sqrt(spdiags(A,0)))',Np,Np);
		inv_sqrt_C   = sparse(1:Nr,1:Nr,(1.0./sqrt(spdiags(C,0)))',Nr,Nr);

		A  = inv_sqrt_A*A  *inv_sqrt_A;
		B1T= inv_sqrt_A*B1T*inv_sqrt_C;
		B2 = inv_sqrt_C*B2 *inv_sqrt_A;
		C  = inv_sqrt_C*C  *inv_sqrt_C;

		W_mat=W_mat*inv_sqrt_C;

		% rhs
		rhsR=zeros(Np+Nr,1);
		rhsR(1:Np)       = inv_sqrt_A * f1;
		rhsR(Np+1:Np+Nr) = ort_proj(inv_sqrt_C*f2,W_mat);
	end


	
	%
	% PREPROCESS TIME
	%
  time_prec=tic;
	time_A=tic;
	
  % define approxiamte inverse of (~A)^{-1}  
  inverseA_approach=controls.inverse11;
	ctrl_inner11 = controls.ctrl_inner11;

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
    		    sprintf('B%sA%d',ctrl_inner11.label,i));
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
	cpu_assembly_inverseA=toc(time_A);
	

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

 
	P = @(rho) ort_proj(rho,W_mat);

	% Define action of preconditoner
	if ( ground_A == 1)
		prec = @(x) SchurCA_based_preconditioner(x, ...
																						 @(y) inv_A(y),...
																						 @(z) P(inverse_block22(P(z))),...
																						 @(y) B1T*P(y),...
																						 @(z) P(B2*z),...
																						 Np,Nr,...
																						 controls.outer_prec,ncellphi);
	else
		prec = @(x) SchurCA_based_preconditioner(x, ...
																						 @(y) ort_proj(inv_A(y),A_kernels),...
																						 @(z) P(inverse_block22(P(z))),...
																						 @(y) B1T*P(y),...
																						 @(z) P(B2*z),...
																						 Np,Nr,...
																						 controls.outer_prec,ncellphi);

	end
	outer_timing=tic;

	[d,info_J]=apply_iterative_solver(@(v) ...
																		 apply_saddle_point(v,@(x) A*x ,...
																												@(y) B1T*y,...
																												@(x) P(B2*x),...
																												@(y) P(C*P(y)),...
																												Np,Nr,N), ...
																			rhsR, controls.ctrl_outer, ...
																			@(z) prec(z),[],controls.left_right );
		
	info_J.print();
	
	if (diagonal_scaling == 1)
		d(1:Np)      = inv_sqrt_A * d(1:Np);
		d(1+Np:Np+Nr)= inv_sqrt_C * d(1+Np:Np+Nr);

		
		sqrt_A   = sparse(1:Np,1:Np,(1./spdiags(inv_sqrt_A,0))',Np,Np);
		sqrt_C   = sparse(1:Nr,1:Nr,(1./spdiags(inv_sqrt_C,0))',Nr,Nr);

		% restore 
		A  = sqrt_A*A  *sqrt_A;
		B1T= sqrt_A*B1T*sqrt_C;
		B2 = sqrt_C*B2 *sqrt_A;
		C  = sqrt_C*C  *sqrt_C;

		W_mat=zeros(N,Nr);
		for k = 1:N
			f_lambda(k) = v' * get_slice(rho,k, ncellrho, N) - 1; % F^k_lambda
			W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = v'; % d F^k_lambda drho_k 
		end	
		
		% rhs
		rhsR(1:Np) = f1;
		rhsR(Np+1:Np+Nr) = ort_proj(f2,W_mat);
	end
										

	d(Np+1:Np+Nr)=ort_proj(d(Np+1:Np+Nr),W_mat);
	
	
	% get s increment
  ds = JF.ss\(-F.s-JF.sr*d(Np+1:Np+Nr));
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

	
	if ( strcmp(inverseA_approach,'full'))
    inner_iter=invA.cumulative_iter;
		inner_nequ=invA.cumulative_nequ;
    inner_cpu=invA.cumulative_cpu;
		invA.kill();
  elseif( strcmp(inverseA_approach,'diag'))
   	inner_nequ=ncellphi;
		inner_iter=0;
		inner_cpu=0;
    for i=1:Nt
      inner_iter=inner_iter+diag_block_invA(i).cumulative_iter;
			inner_cpu=inner_cpu+diag_block_invA(i).cumulative_cpu;
      diag_block_invA(i).kill();
    end
  end

	
	if (strcmp(inverse22,'diagA'))
		inner_nequ2=Nr;
		inner_iter2=inv_SCA.cumulative_iter;
    inv_SCA.kill();
  elseif (strcmp(inverse22,'diagA'))
		inner_nequ2=Nr;
		% set here
		inner_iter2=0;
  end

	resume_msg=sprintf('outer: %d res=%1.1e [%1.1e,%1.1e,%1.1e] iter=%d cpu=%1.1e| A : n=%d it=%d bt=%1.1e  | S : n=%d it=%d bt=%1.1e',...
   									 info_J.flag,resnorm,resp,resr,ress,outer_iter,outer_cpu,...
										 inner_nequ, inner_iter,cpu_assembly_inverseA,...
   									 inner_nequ2, inner_iter2,cpu_assembly_inverseS);

	normd = norm(d);   
  prec_cpu=0;
  inner_nequ3=0;
  inner_iter3=0;

elseif (sol==20)
	index_agmg=0;

	
	[A,B1T,B2,C,f1,f2,W_mat]=get_saddle_point(JF,F);
	B1T_time=-JF.B1T_time;
	B1T_space=-JF.B1T_space;
	P = @(y) P_apply(y,area2h);
	PT = @(y) PT_apply(y,area2h);
	
	rhsR=[f1;P(f2)];
	% assembly SAC=-(C+B2 * diag(A)^{-1} B1T) 
  % reduced to system only in rho
	invDiagA=sparse(1:Np,1:Np,(1./spdiags(A,0))',Np,Np);

	reduced=1;
	primal=1;

	

	ground=controls.ground;
	if (ground)
		% ground system 
		[A,B1T_time,B1T_space,B2_time,B2_space,f1]=ground_saddle(A,B1T_time,B1T_space,f1,N);
		B1T=B1T_time+B1T_space;
		B2 = B2_time+B2_space;
		rhsR=[f1;P(f2)];
	end

	if (controls.diagonal_scaling)
		% scale system by its diagonal
		[sqrt_diagA,sqrt_diagC,inv_sqrt_diagA,inv_sqrt_diagC] = get_diagonal_scaling (A,C);
		
		[A,B1T,B2,C,f1,f2,W_mat]=scale_saddle_point(A,B1T,B2,C,...
																								f1,f2,W_mat,...
																								inv_sqrt_diagA,...
																								inv_sqrt_diagC);
		B1T_time =inv_sqrt_diagA* B1T_time *inv_sqrt_diagC;
		B1T_space=inv_sqrt_diagA* B1T_space *inv_sqrt_diagC;
		B1T=B1T_time+B1T_space;
		invDiagA=sparse(1:Np,1:Np,(1./spdiags(A,0))',Np,Np);

		disp('fix diagonal scaling')
		return
		rhsR=[f1;P(f2)];
	end
	time_A=tic;
	
  % define approxiamte inverse of (~A)^{-1}  
  inverseA_approach=controls.inverse11;
	ctrl_inner11 = controls.ctrl_inner11;
	
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    inverseA=sparse_inverse;
		ctrl_loc=ctrl_solver;
		index_agmg=index_agmg+1;
		ctrl_loc.init(ctrl_inner11.approach,...
    							ctrl_inner11.tolerance,...
    							ctrl_inner11.itermax,...
    							ctrl_inner11.omega,...
    							ctrl_inner11.verbose,...
    							strcat('A',ctrl_inner11.label),...
									index_agmg);
		
    inverseA.init(A+controls.relax4inv11*speye(Np,Np),ctrl_loc);
    inverseA.cumulative_iter=0;
    inverseA.cumulative_cpu=0;

    % define function
    invA = @(x) inverseA.apply(x);

		inner_nequ=Np;

		nsysA=1;
  elseif( strcmp(inverseA_approach,'diag'))
    % partion matrix 
    nAi=ncellphi;
   
    % use block inverse
    diag_block_invA(Nt,1)=sparse_inverse;
    ctrl_loc=ctrl_solver;
   
    for i=1:Nt
			index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
    		    ctrl_inner11.tolerance,...
    		    ctrl_inner11.itermax,...
    		    ctrl_inner11.omega,...
    		    ctrl_inner11.verbose,...
    		    sprintf('%sA%d',ctrl_inner11.label,i),...
						index_agmg);
      diag_block_invA(i).name=sprintf('inverse A%d',i);

      % create local block and passing to solver
      % with a potential relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + ...
	       controls.relax4inv11*speye(nAi,nAi) ;
      diag_block_invA(i).init(matrixAi,ctrl_loc);

			nsysA=Nt;
    end

    % define function
    invA  = @(x) apply_block_diag_inverse(diag_block_invA,x);

		inner_nequ=ncellphi;
  end
	cpu_assembly_inverseA=toc(time_A);


	% define implitely action of the Schur complement
	apply_schur=@(x) P(C*PT(x)+B2*invA(B1T*PT(x)));

	
	%
	% ASSEMBLY SCHUR COMPLEMENT INVERSE
	%
	timeS=tic;
	if (strcmp(controls.inverse22,'diagA'))
		SCA = (B2*invDiagA*B1T);
		if (reduced)
			SCA=SCA+C;
		end
		inv_SCA=sparse_inverse;
		ctrl_inner22 = controls.ctrl_inner22;
		index_agmg=index_agmg+1;
		ctrl_loc=ctrl_solver;
		ctrl_loc.init(ctrl_inner22.approach,...
    							ctrl_inner22.tolerance,...
    							ctrl_inner22.itermax,...
    							ctrl_inner22.omega,...
    							ctrl_inner22.verbose,...
    							'SCA',...
									index_agmg);
		inv_SCA.init(SCA + controls.relax4inv22*speye(Nr,Nr),ctrl_loc);
		inv_SCA.info_inverse.label='schur_ca';
		inv_SCA.cumulative_iter=0;
		inv_SCA.cumulative_cpu=0;


		%define inverse action
		inverseS=	@(z) inv_SCA.apply(z);
		
	elseif (contains(controls.inverse22,'commute'))
		lrb=controls.lrb;
		diagC = spdiags(C,0);
		invC = sparse(1:Nr,1:Nr,(1.0./spdiags(C,0))',Nr,Nr);
		Mphi=JF.Mxt;
		Mrho=-JF.rs;
		inv_Mphi  = sparse(1:Np,1:Np,(1./spdiags(Mphi,0))',Np,Np);
		inv_Mrho  = sparse(1:Nr,1:Nr,(1./spdiags(Mrho,0))',Nr,Nr);
		
		%[pr,A_rho, approx_S,start]=...
		%build_approx_Schur_components( A,B1T,B2,C,invC,lrb,JF,controls.inverse22,controls.mode_Arho);

		inv_Dr = sparse(1:Nr,1:Nr,(1./spdiags(Dr,0))',Nr,Nr);
		
		neit=size(JF.divt_rho,2);
		rho=spdiags(-JF.ss,0);
		rho_edge=spdiags(JF.Rst_rho*rho,0,neit,neit);
		A_rho = - JF.divt_rho*(rho_edge)*JF.gradt_rho;

		
		% rho defined on cells
		RHt = assembleRHt(N,ncellphi);
		vect_rho_on_h = RHt*(JF.It*[F.rho_in;rho;F.rho_f]);
		rho_h =  sparse(1:Np,1:Np,(vect_rho_on_h)',Np,Np);
		inv_rho_h =  sparse(1:Np,1:Np,(1./vect_rho_on_h)',Np,Np);
		

		rec_time = assembleRHt(N-1,ncellrho);
		nei=size(JF.grad,2);
		gradt = assembleGradt(N,ncellphi,nei,JF.div');
		te = size(gradt,1);
		Rst = repmat({JF.Rs},1,N+1);
		Rst = blkdiag(Rst{:});
		I_cellphi_cellrho = assembleIt(N-1,ncellphi,ncellrho,JF.I);
		%  (N+1)*Kh<-(N+1)*e<-(N+1)*e<-(N+1)*Kh <-(N+1)*K2h <- N*K2h
		new_B1T_space = Rst' * spdiags(JF.gradphi,0,te,te) * gradt * I_cellphi_cellrho * rec_time';
		
		
		

		if (0)
			nei=size(JF.gradphi,1)/(N+1);
			for i=1:Nt
				lapl_phi_i=JF.div*JF.gradphi((i-1)*nei+1:i*nei);
				fprintf('%1.1e <= Laplacian phi %d <= %1.1e\n ', min(lapl_phi_i), i, max(lapl_phi_i))
			end
		end

		%S22=B1T_time'*inv_Mphi*B1T_time;
		S22=B2*inv_Mphi*(B1T_time+new_B1T_space);
		
		if (controls.mode_inverse22==1)
			approx_S =  C*inv_Mrho*A_rho + S22;
		elseif (controls.mode_inverse22==2)
			approx_S =  Ds*inv_Mrho*A_rho + Dr*inv_Mrho* S22;
		end
		
		ctrl_inner22 = controls.ctrl_inner22;
		inv_SCA=sparse_inverse;
		ctrl_loc=ctrl_solver;
		index_agmg=index_agmg+1;
		ctrl_loc.init('agmg',...
     							ctrl_inner22.tolerance,...
     							ctrl_inner22.itermax,...
     							0.0,...
     							ctrl_inner22.verbose,...
     							'SCA',...
		 							index_agmg);
		inv_SCA.init(approx_S+1e-10*speye(Nr),ctrl_loc);
		inv_SCA.cumulative_iter=0;
		inv_SCA.cumulative_cpu=0;

		

		%inverseS	= @(y) Dr*inv_Mrho*A_rho*(inv_SCA.apply(y));
		if (controls.mode_inverse22==1)
			inverseS	= @(y) inv_Mrho*A_rho*(inv_SCA.apply(y));
		elseif (controls.mode_inverse22==2)
			inverseS	= @(y) inv_Mrho*A_rho*(inv_SCA.apply(Dr*inv_Mrho*y));
		end
		
	elseif (strcmp(controls.inverse22,'full'))
		schur=prec_times_matrix(apply_schur,speye(Nr,Nr));

		inv_Schur=sparse_inverse;
		ctrl_loc=ctrl_solver;
		index_agmg=index_agmg+1;
		ctrl_loc.init('direct',...
     							1e-12,...
     							200,...
     							0.0,...
     							1,...
     							'SCA',...
		 							index_agmg);
		inv_Schur.init(schur+1e-12*speye(Nr),ctrl_loc);

		inverseS=@(y) inv_Schur.apply(y);

		
	elseif(strcmp(controls.inverse22,'lsc'))
		
		invQ = inv_Mphi;
		%invQ = inv_Mphi*inv_Mphi;
		%invQ = Mphi;
		%invQ = speye(Np);
		% parameters
		gamma=1.0;
		alpha=0e-2;

		
		SCA = (gamma*C+B2*invQ*B1T);
		inv_SCA=sparse_inverse;
		ctrl_inner22 = controls.ctrl_inner22;
		index_agmg=index_agmg+1;
		ctrl_loc=ctrl_solver;
		ctrl_loc.init(ctrl_inner22.approach,...
    							controls.tolerance_preprocess,...
    							ctrl_inner22.itermax,...
    							ctrl_inner22.omega,...
    							ctrl_inner22.verbose,...
    							'SCA',...
									index_agmg);
		inv_SCA.init(SCA + controls.relax4inv22*speye(Nr,Nr),ctrl_loc);
		inv_SCA.info_inverse.label='schur_ca';
		inv_SCA.cumulative_iter=0;
		inv_SCA.cumulative_cpu=0;

		central_block=B2*invQ*A*invQ*B1T;
		SCA = (C+B2*invDiagA*B1T);
		inv_diag_SCA=sparse(1:Nr,1:Nr,(1./spdiags(SCA,0))',Nr,Nr);

		inverse_bbt=@(y) P(inv_SCA.apply(PT(y)));
		
		inverseS=@(y) inverse_bbt(P(central_block*PT(inverse_bbt(y))))+alpha*inv_diag_SCA*y;

	elseif(strcmp(controls.inverse22,'stationary'))
		diagC = spdiags(C,0);
		invC = sparse(1:Nr,1:Nr,(1.0./spdiags(C,0))',Nr,Nr);
		inverseS=@(y) stationary_iterative_methods(apply_schur,y,zeros(Nr,1),...
																							 controls.ctrl_inner22.tolerance,...
																							 controls.ctrl_inner22.itermax,...
																							 controls.ctrl_inner22.verbose,...
																							 @(x) controls.ctrl_inner22.omega*invC*x);
	end
	cpu_assembly_S=toc(timeS);

	%
	% PREPROCESS NULL SPACE METHOD
	%
	inverse_S=@(y) -PT(apply_Schur_inverse(y,inverseS,apply_schur));
	

	if (ground)
		prec = @(x) SchurCA_based_preconditioner(x, @(y) invA(y),...
																						 @(y) inverse_S(y),...
																						 @(y) B1T*y,...
																						 @(z) P(B2*z),...
																						 Np,Nr,...
																						 controls.outer_prec,ncellphi,P,PT);%,@(y) SCA*y);
		
	else
		prec = @(x) SchurCA_based_preconditioner(x, @(y) ort_proj(invA(y),A_kernels),...
																						 @(y) PT(inverse_S(y)),...
																						 @(y) B1T*PT(y),...
																						 @(z) P(B2*z),...
																						 Np,Nr,...
																						 controls.outer_prec,ncellphi,P,PT);%,@(y) SCA*y);
		
	end
	outer_cpu=tic;
	compute_eigen=0;
	if (compute_eigen)
		J=[
			 A, B1T;
			 prec_times_matrix(P,B2), -prec_times_matrix(P,C)
		];
		temp=prec_times_matrix(prec,J);
		%temp=matrix_times_prec(J,prec);

		[real,imag,eigenvec,perm]= study_eigenvalues(temp,'never gonna works');
		
		figure
		plot(abs(real))
		title('prec S')
		%return
	end
		
	if (controls.scaling_rhs)
		scaling=1.0/(norm(rhsR)/norm([f;P(g);h]));
	else
		scaling=1.0;
	end
	
	if (ground)
		rhsR=[f1;P(f2)];
		
		[d,info_J]=apply_iterative_solver(@(v) ...
																			 apply_saddle_point(v,@(x) A*x ,...
																													@(y) B1T*(y),...
																													@(x) P(B2*x),...
																													@(y) P(C*(y)),...
																													Np,Nr,N), ...
																			rhsR, controls.ctrl_outer, ...
																			@(z) prec(z),[],controls.left_right,scaling );
	else
		rhsR=[f1;P(f2)];
		[d,info_J]=apply_iterative_solver(@(v) ...
																			 apply_saddle_point(v,@(x) A*x ,...
																													@(y) B1T*PT(y),...
																													@(x) P(B2*x),...
																													@(y) P(C*PT(y)),...
																													Np,Nr,N), ...
																			rhsR, controls.ctrl_outer, ...
																			@(z) prec(z),[],controls.left_right,scaling);
	end
	
	% get s increment
	d(Np+1:Np+Nr)=PT(d(Np+1:Np+Nr));
	
	if (controls.diagonal_scaling == 1)
		d(1:Np)      = inv_sqrt_diagA * d(1:Np);
		d(1+Np:Np+Nr)= inv_sqrt_diagC * d(1+Np:Np+Nr);
		
		[A,B1T,B2,C,f1,f2,W_mat]=get_saddle_point(JF,F);
		B1T_time=-JF.B1T_time;
		B1T_space=-JF.B1T_space;
		P = @(y) ort_proj(y,W_mat);
		rhsR=[f1;P(f2)];
		disp('fix diagonal scaling')
		return 
	end

	outer_cpu=toc(outer_cpu);
	
	% get ds
	ds = JF.ss\(-F.s-JF.sr*d(Np+1:Np+Nr));
	d = [d; ds];
	
	% get lambda increment dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)
  dl=get_dl(JF,F,d);
  
  % correction of phi increment 
	d(1:Np) = dp_correction(d(1:Np),JF,dl);

	%info_J.print();
	
	%
	% STORE INFO
	%
	outer_iter=uint32(info_J.iter);
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;
  normd = norm(d);

	if ( strcmp(inverseA_approach,'full'))
		inner_iter=inverseA.cumulative_iter;
		inner_cpu=inverseA.cumulative_cpu;
		inverseA.kill();
	elseif ( strcmp(inverseA_approach,'diag'))
		inner_iter=0;
		inner_cpu=0;
		for i=1:Nt
			inner_iter=inner_iter+diag_block_invA(i).cumulative_iter;
			inner_cpu=inner_cpu+diag_block_invA(i).cumulative_cpu;
			diag_block_invA(i).kill();
		end
	end
		
  inner_nequ2=Nr;
	if (  strcmp(controls.inverse22,'diagA') ||...
				contains(controls.inverse22,'commute'))
		inner_iter2=inv_SCA.cumulative_iter;
		inner_nequ2=Nr;
		inv_SCA.kill();

	elseif (  strcmp(controls.inverse22,'lsc') )
		inner_iter2=inv_SCA.cumulative_iter;
		inv_SCA.kill();
	end
	
	cpu_assembly_inverseS=cpu_assembly_S;
	prec_cpu=cpu_assembly_inverseA+cpu_assembly_inverseS;
	% check res
  [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
	relres=resnorm; 
	resume_msg=sprintf('outer: %d res=%1.1e [%1.1e,%1.1e,%1.1e] iter=%d cpu=%1.1e| A %s : in/out=%0.1f bt=%1.1e  | S : in/out=%d bt=%1.1e',...
   									 info_J.flag,resnorm,resp,resr,ress,outer_iter,outer_cpu,...
										 controls.inverse11,(inner_iter/outer_iter)/nsysA,cpu_assembly_inverseA,...
   									 inner_iter2/outer_iter,cpu_assembly_inverseS);


	if (0)
		[A,B1T,B2,C,f1,f2,W_mat]=get_saddle_point(JF,F);
		Sys = [A B1T; B2 -C];
		fp=f1;

		indc=1;
		%[dmax,indc]=max(spdiags(A(1:ncellphi:1:ncellphi),0));
		Sys(indc,:) = sparse(1,Np+Nr);
		Sys(:,indc) = sparse(Np+Nr,1);
		Sys(indc,indc) = 1;
		fp(indc)=0;
		
		d_exact = Sys\[fp;f2];

		d_exact(1:Np)=ort_proj(	d_exact(1:Np),A_kernels);
		d(1:Np)=ort_proj(	d(1:Np),A_kernels);
		d(1+Np:Np+Nr)=ort_proj(	d(1+Np:Np+Nr),W_mat);

		
		figure
		
		plot (d(1:Np+Nr)-d_exact)
		return
	end	
	% grounded_node=400;
	% controls_loc = struct('save_data',0,...
	% 									'indc',grounded_node,...
	% 									'sol',1,...
	% 									'swap_sign',1)
	% approach_string='direct';
	% [omegak, info_solver_newton,norm_ddd,resume_msg] = solvesys(JF,F, controls_loc);


	% figure
	% plot(omegak-d)


	% sww
	% return
elseif sol==220
	% SIMPLE preconditioenr to the full system

	% weights
	Mphi=JF.Mxt;
	Mrho=-JF.rs;
	inv_Mphi  = sparse(1:Np,1:Np,(1./spdiags(Mphi,0))',Np,Np);
	inv_Mrho  = sparse(1:Nr,1:Nr,(1./spdiags(Mrho,0))',Nr,Nr);
	P = @(y) P_apply(y,area2h);
	PT = @(y) PT_apply(y,area2h);
	
  % copy controls
  indc=controls.indc;
  sol=controls.sol;
  ctrl_inner22=controls.ctrl_inner22;
  ctrl_outer=controls.ctrl_outer;
  verbose=controls.verbose;
	time_prec=tic;
	
	% define approxiamte inverse of (~A)^{-1}  
	invdiagA = 1.0./(spdiags(A,0)+controls.relax4inv11);
  invDiagA = sparse(1:Np,1:Np,(invdiagA)',Np,Np);

  % use block inverse
	ctrl_inner11 = controls.ctrl_inner11;
  diag_block_invA(Nt,1)=sparse_inverse;
  ctrl_loc=ctrl_solver;

	nAi=ncellphi;
	time_A=tic;
  for i=1:Nt
    ctrl_loc=ctrl_solver;
		ctrl_loc.init(ctrl_inner11.approach,...
    							ctrl_inner11.tolerance,...
    							ctrl_inner11.itermax,...
    							ctrl_inner11.omega,...
    							ctrl_inner11.verbose,...
    							sprintf('A%d',i));
    diag_block_invA(i).name=sprintf('invA%d',i);

    % create local block and passing to solver
    % with a potential relaxation
    matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + ...
						 controls.relax4inv11*speye(nAi,nAi) ;
    diag_block_invA(i).init(matrixAi,ctrl_loc);     
  end

  % define function
  inv_A  = @(x) apply_block_diag_inverse(diag_block_invA,x);
	cpu_assembly_inverseA=toc(time_A);

	% define implitely action of the Schur complement
	apply_schur=@(x) P(C*PT(x)+B2*inv_A(B1T*PT(x)));




	build_S=tic;
	
	if (contains(controls.approach_Schur_rho,'commute'))
		%
		% assembly Schur_rho=(B2 *inv_Mphi*B1T)*A_rho^{-1}*Mrho
		%                   = S22*A_rho^{-1}*Mrho
		neit=size(JF.divt_rho,2);
		rho=spdiags(-JF.ss,0);
		rho_edge=spdiags(JF.Rst_rho*rho,0,neit,neit);
		A_rho = - JF.divt_rho*(rho_edge)*JF.gradt_rho;
		%A_rho = - JF.divt_rho*JF.gradt_rho;

		% rho defined on cells
		RHt = assembleRHt(N,ncellphi);
		vect_rho_on_h = RHt*(JF.It*[F.rho_in;rho;F.rho_f]);
		rho_h =  sparse(1:Np,1:Np,(vect_rho_on_h)',Np,Np);
		inv_rho_h =  sparse(1:Np,1:Np,(1./vect_rho_on_h)',Np,Np);

		if strcmp( controls.approach_Schur_slack,'dual')
			S22 = (B2*inv_Mphi*B1T);

			% init inverse of full S
			inv_S22=sparse_inverse;
			inv_S22.init( S22 + controls.relax4inv22 * speye(Nr,Nr),ctrl_inner22);
			inv_S22.info_inverse.label='schur22';
			inv_S22.cumulative_iter=0;
			inv_S22.cumulative_cpu=0;
			inverse_block22 = @(x) - inv_Mrho * A_rho * inv_S22.apply(x);

			
			% (S_s)^{-1}= (Dr + Ds  (S22)^{-1}*M)^{-1}
			%           = (Dr M^{-1} S22 S22^{-1}M  + Ds M^{-1} Ar * S22^{-1}*M)^{-1}
			%           = M^{-1} S22 (Dr M^{-1} S22 + Ds M^{-1} Ar)^{-1}
			%           = M^{-1} S22 * S33
			S33  = (Dr*inv_Mrho*S22 + Ds*inv_Mrho*A_rho);

			ctrl_inner33=controls.ctrl_inner33;
			inv_S33 = sparse_inverse;
			inv_S33.init(S33+controls.relax4inv22*speye(Nr,Nr),ctrl_inner33);
			inv_S33.info_inverse.label='schur_ss';
			inv_S33.cumulative_iter=0;
			inv_S33.cumulative_cpu=0;
			inverse_block33 = @(x) inv_Mrho*S22*(inv_S33.apply(x));
			%inverse_block33 = @(x) inv_Mrho*apply_schur(inv_Mrho*S22*(inv_S33.apply(x)));		
			%inverse_block33 = @(y) apply_iterative_solver( @(x) Dr*apply_schur(Mrho*x)+Ds*x, y, ctrl_inner33, @(z) z);
	
		elseif (strcmp( controls.approach_Schur_slack,'primal'))
			if contains(controls.approach_Schur_rho,'Btilde')
				rec_time = assembleRHt(N-1,ncellrho);
				nei=size(JF.grad,2);
				gradt = assembleGradt(N,ncellphi,nei,JF.div');
				te = size(gradt,1);
				Rst = repmat({JF.Rs},1,N+1);
				Rst = blkdiag(Rst{:});
				I_cellphi_cellrho = assembleIt(N-1,ncellphi,ncellrho,JF.I);
				%  (N+1)*Kh<-(N+1)*e<-(N+1)*e<-(N+1)*Kh <-(N+1)*K2h <- N*K2h
				new_B1T_space = Rst' * spdiags(JF.gradphi,0,te,te) * gradt * I_cellphi_cellrho * rec_time';
				
				S22 =  C*inv_Mrho*A_rho + B2*inv_Mphi*(B1T_time+new_B1T_space); %right
			elseif contains(controls.approach_Schur_rho,'Btilde2')
				% define interpolator of rho-functions on phi space
				rec_time = assembleRHt(N,ncellphi);
				It = assembleIt(N,ncell_rho,ncell_phi,JF.I);
				I_all = [sparse(ncell_rho,Nr);speye(Nr,Nr);sparse(ncell_rho,Nr)];
				project_rho2phi = rec_time*It*I_all;

				A_rho = project_rho2phi*A*project_rho2phi;
				S22 =  C*inv_Mrho*A_rho + B2*inv_Mphi*B1T; %right

			elseif contains(controls.approach_Schur_rho,'right')
				S22 =  C*inv_Mrho*A_rho + B2*inv_Mphi*B1T; %right
			elseif contains(controls.approach_Schur_rho,'left')				
				S22 =  A_rho*inv_Mrho*C + B2*inv_Mphi*B1T; %left
			elseif contains(controls.approach_Schur_rho,'both')
				S22 =  A_rho*inv_Mrho*C*inv_Mrho*A_rho + B2*inv_Mphi*A*inv_Mphi*B1T; %both
			end
			%S22 =  C*inv_Mrho*A_rho + B1T_time'*inv_Mphi*B1T_time;
			%S22 =  C*inv_Mrho*A_rho + B1T_time'*inv_Mphi*B1T_time + B1T_space'*inv_Mphi*B1T_space;;
			% init inverse of full S
			inv_S22=sparse_inverse;
			inv_S22.init( S22 + controls.relax4inv22 * speye(Nr,Nr),ctrl_inner22);
			inv_S22.info_inverse.label='schur22';
			inv_S22.cumulative_iter=0;
			inv_S22.cumulative_cpu=0;

			if contains(controls.approach_Schur_rho,'right')
				inverse_block22 = @(x) - inv_Mrho * A_rho * inv_S22.apply(x); %right
			elseif contains(controls.approach_Schur_rho,'Btilde')
				inverse_block22 = @(x) - inv_Mrho * A_rho * inv_S22.apply(x); %right
			elseif contains(controls.approach_Schur_rho,'Btilde2')
				inverse_block22 = @(x) - inv_Mrho * A_rho * inv_S22.apply(x); %right
			elseif contains(controls.approach_Schur_rho,'left')				
				inverse_block22 = @(x) - inv_S22.apply(A_rho *inv_Mrho * x);
			elseif contains(controls.approach_Schur_rho,'both')
				inverse_block22 = @(x) - inv_Mrho * A_rho*inv_S22.apply(A_rho *inv_Mrho * x);  %both
			end
			
			
			ctrl_inner33=controls.ctrl_inner33;
			ctrl_inner33.approach='diag';
			ctrl_inner33.verbose=0;
			inv_S33 = sparse_inverse;
			inv_S33.init(-Dr-controls.relax4inv33*speye(Nr,Nr),ctrl_inner33);
			inv_S33.info_inverse.label='schur_ss';
			inv_S33.cumulative_iter=0;
			inv_S33.cumulative_cpu=0;
			inverse_block33 = @(x) inv_S33.apply(x);
		end
	elseif (strcmp(controls.approach_Schur_rho,'diagA'))
		
		S22 = (B2*invDiagA*B1T);

		inv_S22=sparse_inverse;
		inv_S22.init( S22 + controls.relax4inv22 * speye(Nr,Nr),ctrl_inner22);
		inv_S22.info_inverse.label='schur22';
		inv_S22.cumulative_iter=0;
		inv_S22.cumulative_cpu=0;
		inverse_block22 = @(x) - inv_S22.apply(x);

		% S33^{-1}=(Dr                  +Ds*S^{-1}*M)^{-1}
		%         =(Dr*M^{-1}*S*S^{-1}*M+Ds*S^{-1}*M)^{-1}
		%         =M^{-1} S * (Dr*M^{-1}*S+Ds)^{-1}
		S33  = (Dr*inv_Mrho*S22 + Ds);

		
		inv_S33=sparse_inverse;
		inv_S33.init(S33+controls.relax4inv22*speye(Nr,Nr),ctrl_inner22);
		inv_S33.info_inverse.label='schur_ca';
		inv_S3.cumulative_iter=0;
		inv_S33.cumulative_cpu=0;
		if ( verbose >= 2)
			fprintf('DONE inverse S\n') 
		end 
		inverse_block33 = @(x) inv_Mrho*(S22*(inv_S33.apply(x)));
		%inverse_block33 = @(y) apply_iterative_solver( @(x) Dr*apply_schur(Mrho*x)+Ds*x, y, ctrl_inner33, @(z) z);
	elseif (strcmp(controls.approach_Schur_rho,'lsc'))
		% S^{-1} = (B2 M^-1 B1T)^{-1} B2 M^{-1} A M^{-1} B1T (B2 M^-1 B1T)^{-1}
		S22 = (B2*inv_Mphi*B1T);
		middle_factor=B2*inv_Mphi*A*inv_Mphi*B1T;
		
		inv_S22=sparse_inverse;
		inv_S22.init( S22 + controls.relax4inv22 * speye(Nr,Nr),ctrl_inner22);
		inv_S22.info_inverse.label='Schur_r';
		inv_S22.cumulative_iter=0;
		inv_S22.cumulative_cpu=0;
		inverse_block22 = @(x) -PT(inv_S22.apply(middle_factor*(inv_S22.apply(x))));

		% S33^{-1}=(Dr                  +Ds*S^{-1}*M)^{-1}
		%         =(Dr*M^{-1}*S*S^{-1}*M+Ds*S^{-1}*M)^{-1}
		%         =M^{-1} S * (Dr*M^{-1}*S+Ds)^{-1}
		%S33  = (Dr + Ds);
		inv_Ds  = sparse(1:Nr,1:Nr,(1./spdiags(Ds,0))',Nr,Nr);
		invC = inv_Ds * Dr * inv_Mrho;
		S33 = S22*(invC)*S22+middle_factor;
		
		inv_S33=sparse_inverse;
		inv_S33.init(S33+controls.relax4inv33*speye(Nr,Nr),controls.ctrl_inner33);
		inv_S33.info_inverse.label='Schur_s';
		inv_S3.cumulative_iter=0;
		inv_S33.cumulative_cpu=0;
		if ( verbose >= 2)
			fprintf('DONE inverse S\n') 
		end
		inverse_block33 = @(y) apply_iterative_solver( @(x) Ds*apply_schur(Mrho*x)+Dr*x, y, controls.ctrl_inner33, @(z) z);
		%inverse_block33 = @(x) -inv_Mrho * (inv_S22.apply(inv_S33.apply(inv_S22.apply(inv_Ds*x))));
	end


	cpu_assembly_inverseS=toc(build_S);
	
  % store time required to build prec
  prec_cpu=toc(time_prec);
	
  % set dimension of inner solver
  inner_nequ=Nr;
	inner_nequ2=Nr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % END DUAL SCHUR COMPLEMENT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	
	% Define action of block preconditoner 
	prec11_22 = @(x) SchurCA_based_preconditioner(x,  @(y) ort_proj(inv_A(y),A_kernels),...
																								@(y) PT(inverse_block22(P(y))),...
																								@(y) B1T*PT(y),...
																								@(z) P(B2*z),...
																								Np,Nr,...
																								controls.outer_prec_J,ncellphi);

	
	% global operator
	global_A = @(v) apply_saddle_point(v,@(x) A*x ,...
																		 @(y) B1T*PT(y),...
																		 @(x) P(B2*x),...
																		 @(y) P(-R*PT(y)),...
																		 Np,Nr);
	global_B1T = @(y) [zeros(Np,1);P(M*(y))];
	global_B2  = @(x) Ds*PT(x(Np+1:Np+Nr));
	global_C   = @(z) -Dr*(z);
	
	J=@(x) apply_saddle_point(x,global_A ,...
														global_B1T,...
														global_B2,...
														global_C,...
														Np+Nr,Nr);

	%global rhs
	rhs=[f;P(g);h];
	
	if strcmp( controls.approach_Schur_slack,'dual')
		
		% Define action of block preconditoner 
		prec = @(x) SchurCA_based_preconditioner(x, prec11_22 ,...
																						 @(y) inverse_block33(y),...
																						 global_B1T,...
																						 global_B2,...
																						 Np+Nr,Nr,...
																						 controls.outer_prec,ncellphi);
		
	elseif strcmp( controls.approach_Schur_slack,'primal')
		
		% Define action of block preconditoner 
		prec = @(x) SchurAC_based_preconditioner(x, prec11_22 ,...
																						 @(y) inverse_block33(y),...
																						 global_B1T,...
																						 global_B2,...
																						 Np+Nr,Nr,...
																						 controls.outer_prec,ncellphi);
		
	end

	if (0)
		temp1=prec_times_matrix(J,speye(Np+2*Nr));
		temp=matrix_times_prec(temp1,prec);
		
		figure
		[real,imag,eigenvec,perm]= study_eigenvalues(temp,'never gonna works');
		plot(real,imag,'o')
		
		figure
		semilogy(real)
		title('S')
		return
	end
	
	outer_timing=tic;

	[d,info_J]=apply_iterative_solver(J, rhs, ctrl_outer, @(z) prec(z),[],controls.left_right );
	outer_cpu=toc(outer_timing);



		
	% get lambda increment dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)
  dl=get_dl(JF,F,d);
  
  % correction of phi increment 
	d(1:Np) = dp_correction(d(1:Np),JF,dl);
	
	
	% get info
	outer_iter=uint32(info_J.iter);
	inner_iter=0;
	for i=1:Nt
		inner_iter=inner_iter+diag_block_invA(i).cumulative_iter;
	end
	inner_iter2=inv_S22.cumulative_iter;
	inner_iter3=inv_S33.cumulative_iter;
  
 
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;


 

  relres=res;
  normd = norm(d);
  [ressys,resp,resr,ress]=compute_linear_system_residuum(JF,F,d);


  resume_msg=sprintf('outer: %d ressys=%1.1e [%1.2e,%1.2e,%1.2e] iter=%d cpu=%1.1e| A=%d ; Sr=%d Ss=%d \n',...
   		     info_J.flag,ressys,resp,resr,ress,outer_iter,outer_cpu,...
   		     inner_iter/(outer_iter*Nt),inner_iter2/outer_iter,inner_iter3/outer_iter);

	for i=1:Nt
		diag_block_invA(i).kill();
	end
  inv_S22.kill();
	inv_S33.kill();


elseif (sol==21)
	index_agmg=0;

	
	[A,B1T,B2,C,f1,f2,W_mat]=get_saddle_point(JF,F);
	B1T_time=-JF.B1T_time;
	B1T_space=-JF.B1T_space;

	Q=JF.Mxt;
	inv_Q = sparse(1:Np,1:Np,(1./spdiags(Q,0))',Np,Np);
	
	P = @(y) ort_proj(y,W_mat);
	rhsR=[f1;P(f2)];
	
	ground=0;
	if (ground)
		% ground system 
		[A,B1T_time,B1T_space,B2,f1]=ground_saddle(A,B1T_time,B1T_space,B2,f1,N);
		B1T=B1T_time+B1T_space;		
		rhsR=[f1;P(f2)];
	end


	% define approxiamte inverse of (~A)^{-1}  
	time_A=tic;
  inverseA_approach=controls.inverse11;
	ctrl_inner11 = controls.ctrl_inner11;
	
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    inverseA=sparse_inverse;
		ctrl_loc=ctrl_solver;
		index_agmg=index_agmg+1;
		ctrl_loc.init(ctrl_inner11.approach,...
    							ctrl_inner11.tolerance,...
    							ctrl_inner11.itermax,...
    							ctrl_inner11.omega,...
    							ctrl_inner11.verbose,...
    							strcat('A',ctrl_inner11.label),...
									index_agmg);
		
    inverseA.init(A+controls.relax4inv11*speye(Np,Np),ctrl_loc);
    inverseA.cumulative_iter=0;
    inverseA.cumulative_cpu=0;

    % define function
    invA = @(x) inverseA.apply(x);

		inner_nequ=Np;
		
  elseif( strcmp(inverseA_approach,'diag'))
    % partion matrix 
    nAi=ncellphi;
   
    % use block inverse
    diag_block_invA(Nt,1)=sparse_inverse;
    ctrl_loc=ctrl_solver;
   
    for i=1:Nt
			index_agmg=index_agmg+1;
      ctrl_loc=ctrl_solver;
      ctrl_loc.init(ctrl_inner11.approach,...
    		    ctrl_inner11.tolerance,...
    		    ctrl_inner11.itermax,...
    		    ctrl_inner11.omega,...
    		    ctrl_inner11.verbose,...
    		    sprintf('%sA%d',ctrl_inner11.label,i),...
						index_agmg);
      diag_block_invA(i).name=sprintf('inverse A%d',i);

      % create local block and passing to solver
      % with a potential relaxation
      matrixAi=A((i-1)*nAi+1 :     i*nAi , (i-1)*nAi+1 : i*nAi) + ...
	       controls.relax4inv11*speye(nAi,nAi) ;
      diag_block_invA(i).init(matrixAi,ctrl_loc);     
    end

    % define function
    invA  = @(x) apply_block_diag_inverse(diag_block_invA,x);

		inner_nequ=ncellphi;
  end
	cpu_assembly_inverseA=toc(time_A);
	
	%
	% ASSEMBLY SCHUR COMPLEMENT INVERSE
	%
	timeS=tic;

	BBT=B2*inv_Q*B1T;
	BABT=B2*inv_Q*A*inv_Q*B1T;
	ctrl_inner22 = controls.ctrl_inner22;
	inv_BBT=sparse_inverse;
	ctrl_loc=ctrl_solver;
	index_agmg=index_agmg+1;
	ctrl_loc.init('agmg',...
     						ctrl_inner22.tolerance,...
     						200,...
     						0.0,...
     						ctrl_inner22.verbose,...
     						'SCA',...
		 						index_agmg);

	inv_BBT.init(BBT+1e-10*speye(Nr),ctrl_loc);
	inv_BBT.ctrl.info();
	inv_BBT.cumulative_iter=0;
	inv_BBT.cumulative_cpu=0;

	
	inverse_lsc = @(y) inv_BBT(A*(inv_BBT(y)));


	cpu_assembly_S=toc(timeS);

	%
	% PREPROCESS NULL SPACE METHOD
	%
	cpu_assembly_null_space=0;
	if (controls.null_space)
		time_null=tic;
		invS_Wmat=prec_times_matrix(inverse_lsc,W_mat'); %S^{-1}Z
		MS=W_mat*(invS_Wmat);	 %MS=Z^TA^{-1} Z				
		inverse_MS=sparse_inverse; % create LU factorization of MS
		ctrl_loc=ctrl_solver;
		ctrl_loc.init('direct',...
					 				1e-14,...
					 				200,...
					 				0,...
					 				0,...
					 				'MS',...
					 				84);
		inverse_MS.init(MS,ctrl_loc);
		cpu_assembly_null_space=toc(time_null);

		
		fprintf('preprocess A=%1.1e S=%1.1e NullSpace=%1.1e\n',...
						cpu_assembly_inverseA,cpu_assembly_S,cpu_assembly_null_space) 

		
		%
		% Define action of schur complement
		%
		inv_BBT.ctrl.tolerance=controls.ctrl_inner22.tolerance;	
		null_space_application=@(y) null_space_method(y,P,...
																									inverse_lsc,...
																									W_mat,invS_Wmat,...
																									@(x) inverse_MS.apply(x));
		inverse_S=@(y) -apply_Schur_inverse(y,null_space_application);
	else
		inverse_S=@(y) -apply_Schur_inverse(y,inverse_lsc);
	end


	prec = @(x) SchurCA_based_preconditioner(x, @(y) ort_proj(invA(y),A_kernels),...
																					 inverse_S,...
																					 @(y) B1T*y,...
																					 @(z) P(B2*z),...
																					 Np,Nr,...
																					 controls.outer_prec,ncellphi);

	% get other components changing sign
	M  = -sparse(JF.rs);
	Ds = -sparse(JF.sr);
	Dr = -sparse(JF.ss);
	R  = -sparse(JF.rr);
	R  = -R;

	% get standard component of saddle point
	global_C=-Dr;
			
	global_A = @(x) apply_saddle_point(x,@(y) A*y ,...
																		 @(y) B1T*y,...
																		 @(z) P(B2*z),...
																		 @(z) P(R*z),...
																		 Np,Nr);

	global_B1T = @(y) [zeros(Np,1);P(M*y)];
	global_B2  = @(x) Ds*P(x(Np+1:Np+Nr));
	
	operator=@(d) apply_saddle_point(d,...
																	 global_A ,...
																	 global_B1T,...
																	 global_B2,...
																	 @(y) global_C*y,...
																	 Np+Nr,Nr);
	rhs=[f;P(g);h];


	% SS
	time_ss=tic;


	% invert diagonal compenents
	inv_M   = sparse(1:Nr,1:Nr,(1./spdiags(M,0))',Nr,Nr);			
	inv_Ds = sparse(1:Nr,1:Nr,(1./spdiags(Ds,0))',Nr,Nr);

	% form matrix K=B Ds^{-1} Dr M^{-1} BT + A)
	Schur_s = (B1T*inv_Ds*Dr*inv_M*B2+A);

	% ctrl inverse Schur s
	ctrl_loc=controls.inner33;
	index_agmg=index_agmg+1;
	ctrl_loc.	index_agmg =	index_agmg;
	
	% set inverse
	inv_Schur_s=sparse_inverse;
	inv_Schur_s.init(Schur_s+1e-10*speye(Nr),ctrl_loc);
	inverse_SS=@(z) -(inv_M*B2*inv_Schur_s.apply(B1T*inv_Ds(z)));
	cpu_assembly_ss=toc(time_ss);
	

	% set precondtioenr for full system 
	prec_global = @(x) SchurCA_based_preconditioner(x, prec,...
																									inverse_SS,...
																									global_B1T,...
																									global_B2,...
																									Np+Nr,Nr,...
																									'upper_triang');


	% apply iterative solver
	outer_cpu=tic;
	[d,info_J]=apply_iterative_solver(@(x) operator(x),...
																		rhs, controls.ctrl_outer, ...
																		@(z) prec_global(z),[],controls.left_right );
	
	outer_cpu=toc(outer_cpu);
	


	% get lambda increment dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)
  dl=get_dl(JF,F,d);
  
  % correction of phi increment 
	d(1:Np) = dp_correction(d(1:Np),JF,dl);

	%
	% STORE INFO
	%
	outer_iter=uint32(info_J.iter);
  flag=info_J.flag;
  relres=info_J.res;
  iter=info_J.iter;
  normd = norm(d);

	if ( strcmp(inverseA_approach,'full'))
		inner_iter=inverseA.cumulative_iter;
		inner_cpu=inverseA.cumulative_cpu;
		inverseA.kill();
	elseif ( strcmp(inverseA_approach,'diag'))
		inner_iter=0;
		inner_cpu=0;
		for i=1:Nt
			inner_iter=inner_iter+diag_block_invA(i).cumulative_iter;
			inner_cpu=inner_cpu+diag_block_invA(i).cumulative_cpu;
			diag_block_invA(i).kill();
		end
	end
		
  inner_nequ2=Nr;
	inner_iter2=inv_BBT.cumulative_iter;
	inv_BBT.kill();

	inner_nequ3=Np;
	inner_iter3=inv_Schur_s.cumulative_iter;
	inv_Schur_s.kill();

	
	cpu_assembly_inverseS=cpu_assembly_S+cpu_assembly_null_space;
	prec_cpu=cpu_assembly_inverseA+cpu_assembly_inverseS;
	% check res
  [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d);
	linsol_msg= sprintf('outer: %d res=%1.1e [%1.1e,%1.1e,%1.1e] iter=%d cpu=%1.1e|',...
   									 info_J.flag,resnorm,resp,resr,ress,outer_iter,outer_cpu);
	A_msg= sprintf(' A : n=%d it=%d bt=%1.1e ',...
								 inner_nequ, inner_iter,cpu_assembly_inverseA);

	S_msg= sprintf(' S : n=%d it=%d bt=%1.1e',...
								 inner_nequ2, inner_iter2,cpu_assembly_inverseS);

	SS_msg= sprintf(' SS : n=%d it=%d bt=%1.1e',...
								 inner_nequ3, inner_iter3,cpu_assembly_inverseSS);
	
	resume_msg=strcat(linsol_ms,A_msg,'\n',S_msg,SS_msg);

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
