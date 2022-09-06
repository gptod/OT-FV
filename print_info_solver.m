function [] = print_info_solver(info_solver,fid)
  if (~exist('fid','var') )
    fid=1;
  end
  %fprintf(fid,'flag=%d iter=%d res=%9.4e rhs=%9.4e\n',...
%	  info_solver.flag,info_solver.iter,info_solver.relres,info_solver.rhsnorm);
  fprintf(fid,'outer: %d res %9.4e iter=%d cpu=%9.4e | nequ=%d - inner iter=%d\n',...
	  info_solver.flag,info_solver.relres,info_solver.outer_iter,info_solver.outer_cpu,...
	  info_solver.inner_nequ, info_solver.inner_iter);
  

end
