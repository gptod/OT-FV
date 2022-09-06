function [S_tt,S_tx,S_xx]=assemble_dual_schur_components(A,Nr,N,B1T_time,B1T_space,B1_time,B1_space)
	Np=size(A,1);
	ncellphi=Np/(N+1);
	Nt=N+1;
	
	%
	% define approxiamte inverse of (~A)^{-1}
	%
  inverseA_approach='diag'; % use 'full', or 'diag'
	relax4inv11=1e-12; % system is singular, thus we need to lift it a bit
	% set here other approximate inverse of block11
	ctrl_inner11 = ctrl_solver;
	ctrl_inner11.init('direct',... %approach
										5e-1,... %tolerance
										4,...% itermax
										1e-1,... %omega
										0,... %verbose
										'agmg'); %label
	
  if ( strcmp(inverseA_approach,'full'))
    % set inverse
    inverseA=sparse_inverse;
    inverseA.init(A+relax4inv11*speye(Np,Np),ctrl_inner11);
    inverseA.cumulative_iter=0;
    inverseA.cumulative_cpu=0;

    % define function
    invA = @(x) inverseA.apply(x);
  elseif( strcmp(inverseA_approach,'diag'))
    % partion matrix 
    nAi=ncellphi;
   
    % use block inverse
    diag_block_invA(Nt,1)=sparse_inverse;
    ctrl_loc=ctrl_solver;
   
    for i=1:N+1
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
	       relax4inv11*speye(nAi,nAi) ;
      diag_block_invA(i).init(matrixAi,ctrl_loc);     
    end

    % define functino
    invA  = @(x) apply_block_diag_inverse(diag_block_invA,x);
  end

	
	S_tt=zeros(Nr,Nr);
	S_tx=zeros(Nr,Nr);
	S_xx=zeros(Nr,Nr);

	temp=zeros(Nr,1);

	pb = CmdLineProgressBar('Build Schur components');  
	for i=1:Nr
		pb.print(i,Nr)
		temp(:)=0;
		temp(i)=1;
		
		vt=invA(B1T_time(temp));
		%fprintf('res A=%1.1e \n',norm(A*vt-B1T_time(temp)))
		vx=invA(B1T_space(temp));
		%fprintf('res A=%1.1e \n',norm(A*vx-B1T_space(temp)))
		
		
		S_tt(:,i)=B1_time(vt);
		S_tx(:,i)=B1_time(vx);
		S_xx(:,i)=B1_space(vx);
	end

	S_tt=sparse(S_tt);
	S_tx=sparse(S_tx);
	S_xx=sparse(S_xx);


	if ( strcmp(inverseA_approach,'full'))
		inverseA.kill();
	elseif( strcmp(inverseA_approach,'diag'))
		for i=1:N+1
			diag_block_invA(i).kill();     
		end
	end
