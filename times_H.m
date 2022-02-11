function [B1T_time_out,B1T_space_out,B2_out,C_out,f1_out,f2_out,mat_H] = times_H(B1T_time,B1T_space,B2,C,f1,f2,JF)
	ncellrho = JF.ncellrho;
	N = JF.N;
	
	H=speye(ncellrho-1,ncellrho);   
	H(1:ncellrho-1,ncellrho)=-JF.area2h(1:ncellrho-1)/JF.area2h(ncellrho);
	
	% H=block diagonal (P) 
	mat_H = repmat({H},1,N);
	mat_H = blkdiag(mat_H{:});
	
  B1T_time_out  = B1T_time * mat_H';
	B1T_space_out = B1T_space * mat_H';
  B2_out = mat_H*B2;
	C_out  = mat_H*C*mat_H';

	f1_out=f1;
	f2_out=mat_H*f2;
end

	
