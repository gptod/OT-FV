function [y] =null_space_method(c,P,invS,W_mat,invS_Wmat,inverse_MS)

	%c=f2;
	%invS_Wmat=S\W_mat';
	%MS=W_mat*(invS_Wmat);
	invS_c=invS(c);
	rhs_s=W_mat*invS_c;
	v=inverse_MS(rhs_s);
	%v=MS\rhs_s;
	extra=P(invS_Wmat*v);
	y=P(invS_c)-extra;

	%res=P(S(P(y)))-P(c);
	
	%fprintf('error null space=%1.1e \n ', norm(res)/norm(c))
end
