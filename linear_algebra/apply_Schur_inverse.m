function [y]=apply_Schur_inverse(x,inverse_S,S)


	%fprintf('imbalance g=%1.1e \n ', sum(x))
	y=inverse_S(x);
	
	
	if (exist('S','var') )
		%fprintf('error schur=%1.1e norm(sol)=%1.1e norm(rhs)=%1.1e \n ', ...
		%				norm(S(y)-x)/norm(x),norm(y),norm(x))
	end
