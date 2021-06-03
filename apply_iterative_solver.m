function [sol,info]= apply_iterative_solver(matrix_operator, rhs, ctrl, prec,solini,left_right)

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
  
  if (strcmp(ctrl.approach,'bicgstab'))
    if (strcmp(left_right,'left'))
	      [sol,info.flag,j,info.iter] = bicgstab(@(y) matrix_operator(y),rhs,...
						     ctrl.tolerance,...
						     ctrl.itermax,...
						     @(x) prec(x),[],...% left  prec
						     x0);
    else
      [sol,info.flag,j,info.iter] = bicgstab(@(y) matrix_operator(y),rhs,...
					     ctrl.tolerance,...
					     ctrl.itermax,...
					     [],@(x) prec(x),...% left  prec
					     x0);
    end

	    
	    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
	    %info.flag = (info.res >= ctrl.tolerance);
  elseif (strcmp(ctrl.approach,'pcg'))
      
    [sol,i,j,info.iter] = pcg(@(y) matrix_operator(y),rhs,...
			      ctrl.tolerance,...
			      ctrl.itermax,...
			      @(x) prec(x));%,... 
			      %x0);
    info.flag=i;
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);
	    
  elseif (strcmp(ctrl.approach,'fgmres'))
    [sol,i,res,iter_total] = fgmres(@(y,tol) matrix_operator(y),...
				   rhs,ctrl.tolerance,...
				   'max_iters',ctrl.itermax,...
				   'restart',20,...
				   'x0',x0,...
				   'verb',0,...
				   'P',@(x,tol) prec(x)); %should be right prec

    info.iter=iter_total;
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);

  elseif (strcmp(ctrl.approach,'gmres'))
    if (strcmp(left_right,'left'))
      [sol,i,j,iters] = gmres(@(y) matrix_operator(y),rhs,...
			      20, ... % to fix the total number of iterattion
			      ctrl.tolerance,...
			      int64(ctrl.itermax/20),...
 			      @(x) prec(x),[],... %left ,right
			      x0);      
    else
      [sol,i,j,iters] = gmres(@(y) matrix_operator(y),rhs,...
			      20, ... % to fix the total number of iterattion
			      ctrl.tolerance,...
			      int64(ctrl.itermax/20),...
 			      [],@(x) prec(x),... %left,right
			      x0);
    end
      
    %fprintf('%d %d\n', iters(1),iters(2))
    info.iter=(iters(1)-1)*20+iters(2);
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);

  elseif (strcmp(ctrl.approach,'stationary_iterative_methods'))
    [sol,info.iter,info.res] = stationary_iterative_methods(@(y) matrix_operator(y),rhs,x0,ctrl.tolerance,ctrl.itermax,@(z) prec(z));
    
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);
  else
    disp('IN apply_iterative_solver: solver not defined')
  end
