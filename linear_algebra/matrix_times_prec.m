function [matrix_prec] = matrix_times_prec(matrix,prec)
	N=size(matrix,2);
	matrix_prec=zeros(size(matrix,1),N);
	temp=zeros(N,1);
	for i=1:N
		temp(:)=0;
		temp(i)=1;
		matrix_prec(:,i)=matrix*prec(temp);
	end

end
	
