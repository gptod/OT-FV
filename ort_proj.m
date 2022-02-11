function [x] = ort_proj(x,W_mat)
	% W_mat is assume to be
	% v1^T
	%   v2^T
	%     v3^T
	%       vn^T
	N=size(W_mat,1);
	
	% build projector matrix
	factors=zeros(N,1);
	for k=1:N
		factors(k)=1/norm(W_mat(k,:))^2;
	end
	x=x-W_mat'*(factors.*(W_mat*x));
