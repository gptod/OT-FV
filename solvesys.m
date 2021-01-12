function [d, sol_stat, normd] = solvesys(JF,F,indc,sol)

% solve the linear system d = -JF\F, reducing the system by eliminating the
% additional variable s and imposing the condition d(indc) = 0 to fix the 
% constant of the potential dp


Np = length(F.p);
Nr = length(F.r);
Ns = length(F.s);

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
    
    A = JF.pp; B1 = JF.pr; B2 = JF.rp; invC = -((JF.sr*JF.rs)\JF.ss);
    f1 = F.p;
    f2 = F.r-JF.rs*(JF.ss\F.s);
    S = -A+B1*(invC*B2); S = S(2:end,2:end);
    fp = f1-B1*(invC*f2); fp = fp(2:end);
    dp = S\fp;
    
    sol_stat = struct('flag',0,'relres',[],'iter',[]);
    
    %[dp,flag,res,iter] = agmg(S,fp,0,1e-6,100,0,[],0);
    %sol_stat = struct('flag',flag,'relres',res,'iter',iter);
    
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

end