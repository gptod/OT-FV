function [] = print_info_solver(info_solver)
  fprintf('flag=%d iter=%d res=%9.4e rhs=%9.4e\n',...
	  info_solver.flag,info_solver.iter,info_solver.relres,info_solver.rhsnorm)
end
