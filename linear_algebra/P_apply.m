function [y] = P_apply(x,area2h)
	Nr = size(x);
	ncellrho = size(area2h,1);
	N = Nr / ncellrho;
	inverse_area_domain = 1.0 / sum(area2h);

	y=x;
	for i=1:N
		y((i-1) * ncellrho + 1:(i) * ncellrho) = ...
		y((i-1) * ncellrho + 1:(i) * ncellrho) - ...
		inverse_area_domain *area2h * ...
		sum( 	y((i-1) * ncellrho + 1:(i) * ncellrho) );
	end
end
