function [sol,info]= apply_iterative_solver(matrix_operator, rhs, ctrl, prec, solini,left_right,scaling,prec_right,dim)


	
  if (exist('solini','var') )
    x0=solini;
  else
    x0=zeros(size(rhs,1),1);
  end

  info=info_solver;

  if (~exist('left_right','var') )
    left_right='left';
  else
    if ( strcmp(left_right,'left') || strcmp(left_right,'right') )
      left_right=left_right;
    else
      disp('In apply_iterative_solver string left_right must be left or right')
      return
    end
  end

	if (~exist('scaling','var') )
    scaling=1;
	end

	if (~exist('dim','var') )
    dim=0;
	end

	
  if (strcmp(ctrl.approach,'bicgstab'))
		if (exist('prec_right','var') )
			[sol,infos,j,info.iter] = bicgstab(@(y) matrix_operator(y),rhs,...
																				 ctrl.tolerance,...
																				 ctrl.itermax,...
																				 @(x) prec(x),prec_right,...% left  prec
																				 x0);
		end
    if (strcmp(left_right,'left'))
	      [sol,infos,j,info.iter] = bicgstab(@(y) matrix_operator(y),rhs,...
						     ctrl.tolerance,...
						     ctrl.itermax,...
						     @(x) prec(x),[],...% left  prec
						     x0);
    else
      [sol,infos,j,info.iter] = bicgstab(@(y) matrix_operator(y),rhs,...
					     ctrl.tolerance,...
					     ctrl.itermax,...
					     [],@(x) prec(x),...% right  prec
					     x0);
    end
    info.approach_used='BICGSTAB';
	    
    
  elseif (strcmp(ctrl.approach,'pcg'))
      
    [sol,i,j,info.iter] = pcg(@(y) matrix_operator(y),rhs,...
			      ctrl.tolerance,...
			      ctrl.itermax,...
			      @(x) prec(x),[],... 
			      x0);
    info.flag=i;
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);

    info.approach_used='PCG';
    
  elseif (strcmp(ctrl.approach,'fgmres'))
    nrestart=min(40,ctrl.itermax);
		max_iters=max(int64(ctrl.itermax/nrestart),1);
    [sol,infos,info.res,iter_total] = fgmres(@(y,tol) matrix_operator(y),...
				   rhs,ctrl.tolerance,...
				   'max_iters',max_iters,...
				   'restart',nrestart,...
				   'x0',x0,...
				   'verb',ctrl.verbose,...
					 'tol_exit',ctrl.tolerance,...
				   'P',@(x,tol) prec(x),... %should be right prec
					 'scaling',scaling,...
					 'dim',dim);

    info.iter=iter_total;
    info.flag = (info.res >= ctrl.tolerance);
    info.approach_used='FGMRES';
    
  elseif (strcmp(ctrl.approach,'gmres'))
		nrestart=20;
		disp(scaling)
		correction = 1;%1./scaling;
		if (exist('prec_right','var') && (~isempty( prec_right ))  )
			 [sol,infos,j,iters] = gmres(@(y) matrix_operator(y),rhs,...
			      nrestart, ... % to fix the total number of iterattion
			      ctrl.tolerance*correction,...
			      int64(ctrl.itermax/nrestart),...
 			      @(x) prec(x),prec_right,... %left ,right
			      x0);
		else			
			if (strcmp(left_right,'left'))
				%disp('GMRES LEFT')
				[sol,infos,j,iters] = gmres(@(y) matrix_operator(y),rhs,...
																		nrestart, ... % to fix the total number of iterattion
																		ctrl.tolerance*correction,...
																		int64(ctrl.itermax/nrestart),...
 																		@(x) prec(x),[],... %left ,right
																		x0);      
			else
				%disp('GMRES RIGHT')
				[sol,infos,j,iters] = gmres(@(y) matrix_operator(y),rhs,...
																		nrestart, ... % to fix the total number of iterattion
																		ctrl.tolerance*correction,...
																		int64(ctrl.itermax/nrestart),...
 																		[],@(x) prec(x),... %left,right
																		x0);
			end
		end

    info.approach_used='GMRES';
    %fprintf('%d %d\n', iters(1),iters(2))
    info.iter=(iters(1)-1)*nrestart+iters(2);
    info.flag = (info.res >= ctrl.tolerance);


  elseif ( strcmp(ctrl.approach,'stationary_iterative_methods') )
    if (strcmp(left_right,'left'))
      [sol,info.flag,res,info.iter] = sqmr(@(y) matrix_operator(y),rhs,...
					 ctrl.tolerance,ctrl.itermax,...
					 @(z) prec(z),[],x0)
    else
      [sol,info.flag,res,info.iter] = sqmr(@(y) matrix_operator(y),rhs,...
						  ctrl.tolerance,ctrl.itermax,...
						  [],@(z) prec(z),x0)

    end
    info.approach_used='SQMR'

  elseif (strcmp(ctrl.approach,'stationary_iterative_methods'))
    [sol,info.iter,info.res] = stationary_iterative_methods(@(y) matrix_operator(y),rhs,x0,ctrl.tolerance,ctrl.itermax,@(z) prec(z));
    
    info.flag = (info.res >= ctrl.tolerance);
    info.approach_used='STATIONARY';
  else
    disp('IN apply_iterative_solver: solver not defined')
  end

  ctrl.approach;
  info.flag=infos;
  relres=norm(matrix_operator(sol)-rhs)/norm(rhs);
  info.res=relres;
  info.flag = (info.res >= ctrl.tolerance);
  
  info.rhsnorm=norm(rhs);
  info.balance=sum(rhs);
  info.res=norm(matrix_operator(sol)-rhs);       
  info.realres=info.res/info.rhsnorm;
end
