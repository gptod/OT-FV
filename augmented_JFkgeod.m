function [augJOC] = augmented_JFkgeod (JOC,OC,uk,Number_of_rows, line,integral_constrain,sign,add,verbose) 
%				
% Pass the ORIGINAL (before augmentation) Jacobian and OC 
%
  

  N=JOC.ntimestep;
  ncell=JOC.ncellphi;
  ncell2h=JOC.ncellrho;
		    %				%
  replace=0;

  % store copy of the original
  augJOC=JOC;

  tnr2h = N*ncell2h;
  tnp = (N+1)*ncell;
  phi = uk(1:tnp);
  rho = uk(tnp+1:tnp+tnr2h);
  s = uk(tnp+tnr2h+1:tnp+2*tnr2h);



% k=1,...,N				% 
% G_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s^k)=phi_k^T * Fphi_k(phi_k,        rho_k-1,rho_k            ) +
%                                              rho_k^T * Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k)
% G_{N+1} = deltat * phi^{k+1} I rho^{k}  

%%%%%%%%%%%%%%%%%%%%
% d G / phi
%%%%%%%%%%%%%%%%%%%%
  if (verbose)
     fprintf('ncell=%d \n',ncell)
  end  
  for k=1:Number_of_rows
    % d Gk /d phik
    irow = k ;
    icol = k ;
    vector=get_slice(OC.p,k,ncell,N+1)'+ ...                                                 % Fphi_k
	   get_slice(phi,k,ncell,N+1)' * get_block_jacobian(JOC,'phi','phi',irow,icol) + ... % phi_k*d Fphi_k/d phi_k
	   sign*get_slice(rho,k,ncell2h,N)' * get_block_jacobian(JOC,'rho','phi',irow,icol); % rho_k*d Frho_k/d phi_k
    JOC=set_block_jacobian(JOC,'phi', 'phi', irow, icol,vector, line,add)    
    if (verbose)
      fprintf('norm(Jpp) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
    end
    % d v / d phi_k * Fphi_k + vk * d Fphi_k / Phi_k 
    
    % d Gk /d phi_(k+1)=rho_k^T * d Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k) /d phi_(k+1)
    irow = k;
    icol = k+1;
    vector=sign*get_slice(rho,k,ncell2h,N)' * get_block_jacobian(JOC,'rho','phi',irow,icol) ; %rho_k*d Frho_k/d phi_(k+1)
    JOC=set_block_jacobian(JOC,'phi', 'phi', irow, icol,vector, line,add)
    if (verbose)
      fprintf('norm(Jpp) (%d, %d )=%1.4e\n', irow,icol,norm(vector))
    end
  end

% add additional constraint to remove kernel 
if (integral_constrain)
  % d G_{N+1}/phi_n+1
  irow = N+1;
  icol = N+1;
  augJOC.pp(line+(irow-1)*ncell , 1+(icol-1)*ncell:icol*ncell)=...
  (N+1)* ( Mxt(1:ncell,1:ncell) * It(1:ncell,1:ncell2h))*get_slice(rho,k,ncell2h,N);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d G / rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% G_1(phi_1,phi_2,rho_1,rho_2,s_1)=phi_1^T * Fphi_1(phi_1,      rho_1          ) +
%                                  rho_1^T * Frho_1(phi_1,phi_2,rho_1,rho_2,s_1)
k=1;
nrow=line+(irow-1)*ncell

% d G_1 /d rho1
% irwo icol = lhs position in J.pp
irow = 1 ;
icol = 1 ;
vector=get_slice(phi,k,ncell,N+1)' * get_block_jacobian(JOC,'phi','rho',irow,icol) +... % phi_1 * d Fphi_1/ d rho_1 +
sign*get_slice(OC.r,k,ncell2h,N)'+...                                                 % Frho_1 +
sign*get_slice(rho,k,ncell2h,N)' * get_block_jacobian(JOC,'rho','rho',irow,icol);     % rho_1 * d Frho_1/ d rho_1

JOC=set_block_jacobian(JOC,'phi', 'rho', irow, icol,vector, line,add)
if (verbose)
  fprintf('norm(Jpr) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
end


% d G_1 /d rho{2}
% irwo icol = lhs position in J.pp
k=1;
irow = 1;
icol = 2;
vector=get_slice(phi,k,ncell,N+1)'*get_block_jacobian(JOC,'phi','rho',irow,icol) +...   % phi1 * d Fphi_1/ d rho_2 +
sign*get_slice(rho,k,ncell2h,N)'*get_block_jacobian(JOC,'rho','rho',irow,icol);       % rho1 * d Frho_1/ d rho_2

JOC=set_block_jacobian(JOC,'phi', 'rho', irow, icol,vector, line,add)
if (verbose)
  fprintf('norm(Jpr) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
end



for k=2:min(Number_of_rows,N-1)
% G_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s^k)=phi_k^T * Fphi_k(phi_k,        rho_k-1,rho_k            ) +
%                                              rho_k^T * Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k)

  % d Gphik /d rho{k-1} = d(rho_k^T * Frho_k)
  % lower diagonal
  irow = k ;
  icol = k-1;
  vector=...
  get_slice(phi,k,ncell,N+1)'      * get_block_jacobian(JOC,'phi','rho',irow,icol) +...% phi_k * d Fphi_k/ d rho_{k-1} +
  sign*get_slice(rho,k,ncell2h,N)' * get_block_jacobian(JOC,'rho','rho',irow,icol);     % rho_k * d Frho_k/ d rho_{k-1}
  
  JOC=set_block_jacobian(JOC,'phi', 'rho', irow, icol,vector, line,add)


  if (verbose)
    fprintf('norm(Jpr) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
  end

  % d Gphik /d rhok
  % diagonal
  irow = k ;
  icol = k ;
  vector=...
  get_slice(phi,k,ncell,N+1)' * get_block_jacobian(JOC,'phi','rho',irow,icol) +...   % phi_k * d Fphi_k/ d rhok +
  sign*get_slice(OC.r,k,ncell2h,N)'+...                                                       % Frho_k +
  sign*get_slice(rho,k,ncell2h,N)' * get_block_jacobian(JOC,'rho','rho',irow,icol);       % rho_k * d Frho_k/ d rho_k

  JOC=set_block_jacobian(JOC,'phi', 'rho', irow, icol,vector, line,add)
								      
  if (verbose)
    fprintf('norm(Jpr) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
  end
  
  % d Gphik /d rho{k+1}
  % irwo icol = lhs position in J.pp
  irow = k ;
  icol = k+1;
  vector=...
  get_slice(phi,k,ncell,N+1)'*get_block_jacobian(JOC,'phi','rho',irow,icol) + ... % phi_k* d Fphi_k/ d rho_{k+1}=0 +
  sign*get_slice(rho,k,ncell2h,N)'*get_block_jacobian(JOC,'rho','rho',irow,icol);      % rho_k* d Frho_k/ d rho_{k+1}

  JOC=set_block_jacobian(JOC,'phi', 'rho', irow, icol,vector, line,add)

  if (verbose)
    fprintf('norm(Jpr) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
  end
end

if (Number_of_rows == N)
  k=N;
  
  % d Gphik /d rho{k-1}
  % irwo icol = lhs position in J.pp
  irow = k ;
  icol = k-1;
  vector=get_slice(phi,k,ncell,N+1)'*get_block_jacobian(JOC,'phi','rho',irow,icol) +...  % phi_k * d Fphi_k/ d rho_{k-1} +
  sign*get_slice(rho,k,ncell2h,N)'*get_block_jacobian(JOC,'rho','rho',irow,icol);       % rho_k * d Frho_k/ d rho_{k-1}

  JOC=set_block_jacobian(JOC,'phi', 'rho', irow, icol,vector, line,add)
  if (verbose)
    fprintf('norm(Jpr) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
  end
  
  % d Gphik /d rhok
  % irwo icol = lhs position in J.pp
  irow = k; 
  icol = k;
  vector=get_slice(phi,k,ncell,N+1)' * get_block_jacobian(JOC,'phi','rho',irow,icol) +... % phi_k * d Fphi_k/ d rhok +
	 sign*get_slice(OC.r,k,ncell2h,N)'+...                                                     % Frho_k +
	 sign*get_slice(rho,k,ncell2h,N)' * get_block_jacobian(JOC,'rho','rho',irow,icol);     % rho_k * d Frho_k/ d rho_k

  JOC=set_block_jacobian(JOC,'phi', 'rho', irow, icol,vector, line,add)
  if (verbose)
    fprintf('norm(Jpr) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
  end
end

if (integral_constrain)
  % d GN+1 / d rho_k
  irow = N+1;
  icol = N;
  augJOC.pr(line+(irow-1)*ncell,1+(icol-1)*ncell2h:icol*ncell2h)=...
  (N+1)*get_slice(phi,N+1,ncell,N+1)'* ( Mxt(1:ncell,1:ncell) * It(1:ncell,1:ncell2h));
end

%%%%%%%%%%%%%
% d G / s
%%%%%%%%%%%%
for k=1:Number_of_rows
% G_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s^k)=phi_k^T * Fphi_k(phi_k,        rho_k-1,rho_k            ) +
%                                              rho_k^T * Frho_k(phi_k,phi_k+1,rho_k-1,rho_k,rho_k+1,s_k)
% thus 
% dG/sk =  rho_k^T *  d(Frho_k) /ds^k
  irow = k;
  icol = k;
  vector=sign*get_slice(rho,k,ncell2h,N)' * get_block_jacobian(JOC,'rho','s',irow,icol);
  JOC=set_block_jacobian(JOC,'phi', 's', irow, icol,vector, line,add)

  if (verbose)
    fprintf('norm(Jps) (%d, %d )=%1.2e\n', irow,icol,norm(vector))
  end
end

