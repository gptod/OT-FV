function [x] = zero_mean(x,area)
	% W_mat is assume to be
	% v1^T
	%   v2^T
	%     v3^T
	%       vn^T
	ncell = size(area,1);
	N = size(x,1)/ ncell;
	area_domain = sum(area);
	
	% build projector matrix
	for k = 1 : N
		imbalance = x(1+(k-1)*ncell : k*ncell )'*area/area_domain
		x(1+(k-1)*ncell : k*ncell ) = x(1+(k-1)*ncell : k*ncell ) - imbalance; 
	end
