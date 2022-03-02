function [y] =null_space_method(c,P,PT,invS,W_mat,invS_Wmat,inverse_MS)

	invS_c=invS(c);
	rhs_s=W_mat*invS_c;
	% W^T S^{-1} W v = S^{-1} c
	v=inverse_MS(rhs_s);
	
	% 
	extra=invS_Wmat*v;
	
	y=PT(P(invS_c-extra));
	
	%fprintf('error null space=%1.1e \n ', norm(res)/norm(c))
end
