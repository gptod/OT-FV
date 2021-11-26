function [JOC] = JFkgeod(N,Dt,divt,Mxt,Mxt2h,gradt,It,rhos,drhos,ddrhosa,uk,OC,augmented)

% Compute jacobian of the system of first-order optimality conditions

Nt = N+1;
ncell2h = size(Mxt2h,1)/Nt;
ncell = size(Mxt,1)/Nt;
nei = size(divt,2)/Nt;
tnr2h = N*ncell2h;
tnp = Nt*ncell;
te = Nt*nei;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);
s = uk(tnp+tnr2h+1:tnp+2*tnr2h);

I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
gradphi = gradt*phi;

% derivatives of the optimality conditions in phi
JOCpp = divt*spdiags(rhos,0,te,te)*gradt;
JOCpr = Mxt*Nt*Dt*It*I_all+divt*spdiags(gradphi,0,te,te)*drhos;
JOCps = sparse(tnp,tnr2h);

% derivatives of the optimality conditions in rho
%JOCrp = (Nt*Mxt*Dt*It*I_all)'-drhos'*Mst*(spdiags(gradphi,0,te,te)*gradt);
JOCrr = -0.5*ddrhosa;
JOCrs = -Mxt2h(1:tnr2h,1:tnr2h);

% optimality conditions in s
JOCsp = sparse(tnr2h,tnp);
JOCsr = -spdiags(s,0,tnr2h,tnr2h);
JOCss = -spdiags(rho,0,tnr2h,tnr2h);




JOC.pp = JOCpp; JOC.pr = JOCpr; JOC.ps = JOCps;
JOC.rp = JOCpr'; JOC.rr = JOCrr; JOC.rs = JOCrs;
JOC.sp = JOCsp; JOC.sr = JOCsr; JOC.ss = JOCss;

JOC.ntimestep=N;
JOC.ncellrho = ncell2h;
JOC.ncellphi = ncell;

JOC.DtT=Nt*Dt*It*I_all;
JOC.B1T_time=Mxt*Nt*Dt*It*I_all;
JOC.B1T_space=divt*spdiags(gradphi,0,te,te)*drhos;
  
JOC.LT=spdiags(gradphi,0,te,te)*drhos;
JOC.diagrho=spdiags(rhos,0,te,te);



if ~exist('augmented', 'var') || isempty(augmented)
  augmented=0;
end

if (augmented)		    %				%
  replace=0;
				% 
  % G_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s^k)=phi_k^T * Fphi_k(phi_k,        rho_k-1,rho_k            ) +
  %                                              rho_k^T * Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k)
  % k=1,...,N
  % G_{N+1} = deltat * phi^{k+1} I rho^{k}  

  %%%%%%%%%%%%%%%%%%%%
  % d G / phi
  %%%%%%%%%%%%%%%%%%%%
  for i=1:N
    % copy vectors
    phik  = phi (1+(i-1)*ncell  :(i)*ncell  );
    rhok  = rho (1+(i-1)*ncell2h:(i)*ncell2h);

    Fphik = OC.p(1+(i-1)*ncell  :(i)*ncell  );
    Frhok = OC.r(1+(i-1)*ncell2h:(i)*ncell2h);

    % d Gk /d phik
    % irwo icol = lhs position in J.pp
    irow = i ;
    icol = i ;

    % put in the last row
    JOC.pp(ncell+(irow-1)*ncell , 1+(icol-1)*ncell:icol*ncell) = ...
    Fphik'+ ...                                                 % Fphi_k
    phik' * get_block_jacobian(JOC,'phi','phi',irow,icol) + ... % phi_k*d Fphi_k/d phi_k
    rhok' * get_block_jacobian(JOC,'rho','phi',irow,icol);      % rho_k*d Frho_k/d phi_k

    % d Gk /d phi_(k+1)=rho_k^T * d Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k) /d phi_(k+1)
    irow = i;
    icol = i+1;
    JOC.pp(ncell+(irow-1)*ncell , 1+(icol-1)*ncell:icol*ncell) = ...
    rhok' * get_block_jacobian(JOC,'rho','phi',irow,icol) ; %rho_k*d Frho_k/d phi_(k+1)
    

  end

  % d G_{N+1}/phi_n+1
  rhok  = rho (1+(N-1)*ncell2h:(N)*ncell2h);

  irow = N+1;
  icol = N+1;
  JOC.pp(ncell+(irow-1)*ncell , 1+(icol-1)*ncell:icol*ncell)=...
  (N+1)* ( Mxt(1:ncell,1:ncell) * It(1:ncell,1:ncell2h))*rhok;


	   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % d G / rho
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % G_1(phi_1,phi_2,rho_1,rho_2,s_1)=phi_1^T * Fphi_1(phi_1,rho_1) +
  %                                            rho_1^T * Frho_1(phi_1,phi_2,rho_1,rho_2,s_1)
  i=1;
      % copy vectors
  phik  = phi (1+(i-1)*ncell  :(i)*ncell  );
  rhok  = rho (1+(i-1)*ncell2h:(i)*ncell2h);
  
  Fphik = OC.p(1+(i-1)*ncell  :(i)*ncell  );
  Frhok = OC.r(1+(i-1)*ncell2h:(i)*ncell2h);
  
  % d Gphik /d rho1
  % irwo icol = lhs position in J.pp
  irow = 1 ;
  icol = 1 ;
  JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
  phik' * get_block_jacobian(JOC,'phi','rho',irow,icol) +...   % phi_1 * d Fphi_1/ d rho_1 +
  Frhok'+...                                                % Frho_1 +
  rhok' * get_block_jacobian(JOC,'rho','rho',irow,icol);     % rho_1 * d Frho_1/ d rho_1

  % d Gphik /d rho{2}
  % irwo icol = lhs position in J.pp
  irow = 1;
  icol = 2;
  JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
  phik'*get_block_jacobian(JOC,'phi','rho',irow,icol) +...   % phi1 * d Fphi_1/ d rho_2 +
  rhok'*get_block_jacobian(JOC,'rho','rho',irow,icol);       % rho1 * d Frho_1/ d rho_2

  for i=2:N-1
    % G_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s^k)=phi_k^T * Fphi_k(phi_k,        rho_k-1,rho_k            ) +
    %                                              rho_k^T * Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k)
    % copy vectors
    phik  = phi (1+(i-1)*ncell  :(i)*ncell  );
    rhok  = rho (1+(i-1)*ncell2h:(i)*ncell2h);

    Fphik = OC.p(1+(i-1)*ncell  :(i)*ncell  );
    Frhok = OC.r(1+(i-1)*ncell2h:(i)*ncell2h);
    
    % d Gphik /d rho{k-1} = d(rho_k^T * Frho_k)
    % irwo icol = lhs position in J.pp
    irow = i ;
    icol = i-1;
    JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
    phik'*get_block_jacobian(JOC,'phi','rho',irow,icol) +...   % phi_k * d Fphi_k/ d rho_{k-1} +
    rhok'*get_block_jacobian(JOC,'rho','rho',irow,icol);        % rho_k * d Frho_k/ d rho_{k-1}
    
  
    % d Gphik /d rhok
    % irwo icol = lhs position in J.pp
    irow = i ;
    icol = i ;
    JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
    phik' * get_block_jacobian(JOC,'phi','rho',irow,icol) +...   % phi_k * d Fphi_k/ d rhok +
    Frhok'+...                                               % Frho_k +
    rhok' * get_block_jacobian(JOC,'rho','rho',irow,icol);        % rho_k * d Frho_k/ d rho_k

    % d Gphik /d rho{k+1}
    % irwo icol = lhs position in J.pp
    irow = i ;
    icol = i+1;
    JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
    phik'*get_block_jacobian(JOC,'phi','rho',irow,icol) + ...   % phi_k * d Fphi_k/ d rho_{k+1} +
    rhok'*get_block_jacobian(JOC,'rho','rho',irow,icol);       % rho_k * d Frho_k/ d rho_{k+1}
  end

  i=N;

  % copy vectors
  phik  = phi (1+(i-1)*ncell  :(i)*ncell  );
  rhok  = rho (1+(i-1)*ncell2h:(i)*ncell2h);
  
  Fphik = OC.p(1+(i-1)*ncell  :(i)*ncell  );
  Frhok = OC.r(1+(i-1)*ncell2h:(i)*ncell2h);

  % d Gphik /d rho{k-1}
  % irwo icol = lhs position in J.pp
  irow = i ;
  icol = i-1;
  JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
  phik'*get_block_jacobian(JOC,'phi','rho',irow,icol) +...  % phi_k * d Fphi_k/ d rho_{k-1} +
  rhok'*get_block_jacobian(JOC,'rho','rho',irow,icol);       % rho_k * d Frho_k/ d rho_{k-1}

  % d Gphik /d rhok
  % irwo icol = lhs position in J.pp
  irow = i; 
  icol = i;
  JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
  phik' * get_block_jacobian(JOC,'phi','rho',irow,icol) +...  % phi_k * d Fphi_k/ d rhok +
  Frhok'+...                                                                                % Frho_k +
  rhok' * get_block_jacobian(JOC,'rho','rho',irow,icol);     % rho_k * d Frho_k/ d rho_k

  % d GN+1 / d rho_k
  phik=phi(1+N*ncell:(N+1)*ncell);
  
  irow = N+1;
  icol = N;
  JOC.pr(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
  (N+1)*phik'* ( Mxt(1:ncell,1:ncell) * It(1:ncell,1:ncell2h));

  % d G / s
  for i=1:N
    %
    % G_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s^k)=phi_k^T * Fphi_k(phi_k,        rho_k-1,rho_k            ) +
    %                                              rho_k^T * Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k)
    % thus 
    % dG/sk =  rho_k^T *  d(Frho_k) /ds^k

    % copy vectors
    phik  = phi (1+(i-1)*ncell  :(i)  *ncell  );
    rhok  = rho (1+(i-1)*ncell2h:(i)  *ncell2h);

    Fphik = OC.p(1+(i-1)*ncell  :(i-1)*ncell  );
    Frhok = OC.r(1+(i-1)*ncell2h:(i-1)*ncell2h);

    
    irow = i;
    icol = i;
    JOC.ps(ncell+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
    rhok' * get_block_jacobian(JOC,'rho','s',irow,icol);   
  end
  
end
  

JOC.It=It;

JOC.Mxt=Mxt;


end
