    %B2_invA_B1T=dense_S-C;

    rho=-spdiags(JF.ss,0);
    fprintf('%1.4e<= rho<=%1.4e\n', min(rho),max(rho))


    %B2_invA_B1T=B2*(A\B1T);
    only_time=zeros(Nr,Nr);
    only_space=zeros(Nr,Nr);
    perturbation=zeros(Nr,Nr);
    %global B2_invA_B1T,only_time, only_space,perturbation;
   

    if (0)
      fprintf('%4.2f\n',0)
      for i=1:N
	for j=1:ncellphi		
	  k=(i-1)*ncellphi+j;
	  if (mod(k,100)==0)
	    fprintf('\b\b\b\b\b\b%4.2f ',k/Nr*100)	  
	  end
	  vector=(full(B2_time(k,:)))';
				%size(vector)
				%size(	real_inverseA.matrix)
	  vect_out=real_inverseA.apply(vector);     
	  only_time(:,k)=B2_time*vect_out;
	end
      end
    else
      invA_B2_time=full(A\(B2_time'));
    end
    
    if (0)
      fprintf('%4.2f\n',0)
      for i=1:N
	for j=1:ncellphi		
	  k=(i-1)*ncellphi+j;
	  if (mod(k,100)==0)
	    fprintf('\b\b\b\b\b\b%4.2f ',k/Nr*100)	  
	  end
	  vector=full(B1T_perturbation(:,k));
				%size(vector)
				%size(	real_inverseA.matrix)
	  vect_out=real_inverseA.apply(vector);
	  perturbation(:,k)=B2*vect_out;
	end
      end
    else
      invA_P=full(A\(B1T_perturbation));
    end

    invA_B2_space=full(A\(B2_space)');

    %M=B2_time+B2_space;
    %Mcut=M(1:Nr,ncellpho+1:Np);
    
    %approx=full(Mcut*(A(ncellphi+1:Np,ncellphi+1:Np)\(Mcut')))

    % block_diag_S=zeros(Nr,Nr);
    % for i=1:N
    %   block_diag_S((i-1)*ncellrho+1:i*ncellrho,(i-1)*ncellrho+1:i*ncellrho)=...
    %   dense_S((i-1)*ncellrho+1:i*ncellrho,(i-1)*ncellrho+1:i*ncellrho);
    % end

    % block_tridiag_S=zeros(Nr,Nr);
    % for i=1:N
    %   for j=1:N
    % 	block_tridiag_S((j-1)*ncellrho+1:j*ncellrho,(i-1)*ncellrho+1:i*ncellrho)=...
    % 	diag(diag(dense_S((j-1)*ncellrho+1:j*ncellrho,(i-1)*ncellrho+1:i*ncellrho)));
    %   end
    % end

    % size(JF.B_space)
    % fprintf('%4.2f\n',0)
    % for i=1:N
    %   for j=1:ncellphi		
    % 	k=(i-1)*ncellphi+j;
    % 	if (mod(k,100)==0)
    % 	  fprintf('\b\b\b\b\b\b%4.2f ',k/Nr*100)	  
    % 	end
    % 	vector=(full(JF.B_space(k,:)))';
    % 	vect_out=real_inverseA.apply(vector);
    % 	only_space(:,k)=JF.B_space*vect_out;
    %   end
    % end
    
    %only_space=B2_invA_B1T-only_time-perturbation;
    ctrl_loc=ctrl_solver;
    ctrl_loc.init('direct',1e-14,0,0.0,0,'invS');
    real_inverseS=sparse_inverse;
    real_inverseS.init(dense_S,ctrl_loc);

    M_C            = prec_times_matrix(@(x) real_inverseS.apply(x), C);
    M_time_time   = prec_times_matrix(@(x) real_inverseS.apply(x), B2_time*invA_B2_time);
    M_space_space = prec_times_matrix(@(x) real_inverseS.apply(x), B2_space*invA_B2_space);
    M_time_space  = prec_times_matrix(@(x) real_inverseS.apply(x), B2_time*invA_B2_space);
    M_space_time  = prec_times_matrix(@(x) real_inverseS.apply(x), B2_space*invA_B2_time);
    M_P           = prec_times_matrix(@(x) real_inverseS.apply(x), B2*invA_P);

    

    mc=norm(M_C)
    mtt=norm(M_time_time)
    mts=norm(M_time_space)
    mst=norm(M_space_time)
    mss=norm(M_space_space)
    mp=norm(M_P)
    candidate=M_C+M_time_time+M_space_time;
    norm(candidate-speye(Nr,Nr))


    eigenvalues=study_eigenvalues(candidate, 'S^{-1}( C +B2*invA*B_time)',1);

    figure
    plot(eigenvalues,'o')


    
    tot=M_C+M_time_time+M_space_space+M_time_space+M_space_time+M_P-speye(Nr,Nr);
    norm(tot)
	     


    %spy(tot)
    
    return
    
    % eigenvalues=eig(only_space);
    % real_part=real(eigenvalues);
    % imag_part=imag(eigenvalues);
    % abs_eig=abs(eigenvalues);
    % fprintf('real part eigen (S_space) : min=%1.4e max=%1.4e\n', min(real_part),max(real_part))
    % fprintf('imag part eigen (S_space) : min=%1.4e max=%1.4e\n', min(imag_part),max(imag_part))
    % fprintf('abs       eigen (S_space) : min=%1.4e max=%1.4e\n', min(abs_eig),max(abs_eig))

    % study_eigenvalues(dense_S\(C+B2_time*invA_B2_time)-speye(Nr,Nr), 'S^{-1}( C +B_time * invA * B_timeT)-I',controls.logID);
    % study_eigenvalues(dense_S\(C+B2_time*invA_B2_space)-speye(Nr,Nr), 'S^{-1} (C +B_time * invA * B_space)-I',controls.logID);
    % study_eigenvalues(dense_S\(C+B2_time*invA_P), 'S^{-1} (C+B_time * invA * P)',controls.logID);
    % study_eigenvalues(dense_S\(C+B2_space*invA_B2_time)-speye(Nr,Nr), 'S^{-1} (C+B_space * invA * B_time)-I',controls.logID);
    %study_eigenvalues(dense_S\(B2_space*invA_B2_space+B2_time*invA_B2_space), 'S^{-1} (C+B_space * invA * B_space)-I',controls.logID);
    % study_eigenvalues(dense_S\(C+B2_space*invA_P), 'S^{-1} (C+B_space * invA * P)',controls.logID);
    % study_eigenvalues(dense_S\(C+(B2_time+B2_space) * ( invA_B2_time + invA_P + invA_B2_space) ), 'check',controls.logID);
    %study_eigenvalues((dense_S\(C+B2_time*invA_B2_time+B2_space*invA_B2_space)), 'S^{-1}( C +B_time * invA * B_timeT+B_space * invA * B_spaceT)',controls.logID);
    %saveas(gcf,strcat(controls.basename,'invStimespace.png'));

    %study_eigenvalues((dense_S\(block_diag_S)), 'S^{-1}(blockdiag )',controls.logID);

	    %    saveas(gcf,strcat(controls.basename,'invSdiag.png'));


    %study_eigenvalues((dense_S\(block_tridiag_S)), 'S^{-1}(blocktridiag )',controls.logID);

    %saveas(gcf,strcat(controls.basename,'invStridiag.png'));
  

    d=zeros(Np+2*Nr);
    
    sol_stat = struct('flag',0,...
		      'relres',0,...
		      'ressol',0,...
		      'iter',0,...
		      'rhsnorm',norm([F.p;F.r;F.s]),...
		      'inner_nequ', 0,...
		      'inner_iter', 0,...
		      'outer_iter', 0,...
		      'outer_cpu',  0 );
    return
    
