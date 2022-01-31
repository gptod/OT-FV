function [A,B1T_time,B1T_space,B2,C,f1,f2,H_mat] = times_H(A,B1T_time,B1T_space,B2,C,f1,f2,JF)
	ncellrho = JF.ncellrho;
	N = JF.N;
	
	H=speye(ncellrho-1,ncellrho);   
	H(1:ncellrho-1,ncellrho)=-JF.area2h(1:ncellrho-1)/JF.area2h(ncellrho);
	
	% H=block diagonal (P) 
	mat_H = repmat({H},1,N);
	mat_H = blkdiag(mat_H{:});
	
	A=A;
  B1T_time  = B1T_time * mat_H';
	B1T_space = B1T_space * mat_H';
  B2 = mat_H*B2;
	C  = mat_H*C*mat_H';

	f1=f1;
	f2=mat_H*f2;
end

	
