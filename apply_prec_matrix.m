function [prec_matrix] = apply_prec_matrix(prec,matrix)
	N=size(matrix,1);
	prec_matrix=zeros(N,N);
	for i=1:N
		prec_matrix(:,i)=prec(matrix(:,i));
	end

end
	
