function [sol,info]= apply_iterative_solver(matrix_operator, rhs, ctrl, prec,solini)

  if (exist('solini','var') )
    x0=solini;
  else
    x0=zeros(size(rhs,1),1);
  end

  info=info_solver;
  
  if (strcmp(ctrl.approach,'bicgstab'))
	    [sol,i,j,info.iter] = bicgstab(@(y) matrix_operator(y),rhs,...
					     ctrl.tolerance,...
					     ctrl.itermax,...
					     @(x) prec(x));
	    
	    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
	    info.flag = (info.res >= ctrl.tolerance);
	    
  elseif (strcmp(ctrl.approach,'fgmres'))
    [sol,i,res,iter_total] = fgmres(@(y,tol) matrix_operator(y),...
				   rhs,ctrl.tolerance,...
				   'max_iters',ctrl.itermax,...
				   'restart',20,...
				   'verb',0,...
				   'P',@(x,tol) prec(x));

    info.iter=iter_total;
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);

  elseif (strcmp(ctrl.approach,'gmres'))
    [sol,i,j,iters] = gmres(@(y) matrix_operator(y),rhs,...
			    20, ctrl.tolerance,...
			    ctrl.itermax,...
 			    @(x) prec(x));
    info.iter=iters(1);
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);

  elseif (strcmp(ctrl.approach,'stationary_iterative_methods'))
    [sol,info.iter,info.res] = stationary_iterative_methods(@(y) matrix_operator(y),rhs,x0,ctrl.tolerance,ctrl.itermax,@(z) prec(z));
    
    info.res=norm(matrix_operator(sol)-rhs)/norm(rhs);
    info.flag = (info.res >= ctrl.tolerance);
  else
    disp('IN apply_iterative_solver: solver not defined')
  end
