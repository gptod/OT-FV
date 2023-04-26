function [prec] = assemble_prec_dual(A,B1T,B2,C,tol)
	% assembly SAC=-(C+B2 * diag(A)^{-1} B1T)
	Np=size(A,1);
	Nr=size(C,1);
	
	invDiagA   = sparse(1:Np,1:Np,(1.0./spdiags(A,0))',Np,Np);
	approx_SCA = (C+B2*invDiagA*B1T);

	ctrlS=ctrl_solver;
	ctrlS.init('agmg',...
						 tol,...
						 1000,...
						 0,...
						 1,...
						 'agmg',3);

	inv_SCA=sparse_inverse;
	inv_SCA.init(approx_SCA + 0*speye(Nr,Nr),ctrlS);
	inv_SCA.info_inverse.label='schur_ca';
	prec = @(y) inv_SCA.apply(y);
end
