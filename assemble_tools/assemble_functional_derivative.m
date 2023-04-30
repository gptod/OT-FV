function [out]=assemble_functional_derivative(fun_dfun_ddfun,order,rho,rho_initial,rho_final,area2h)
	% return the value, gradient and hessian of the a rho functional added to the kinetic energy
	% order = 0 : out scalar
	% order = 1 : out vector
	% order = 2 : out matrix

	% type = 'otp' no extra term
	Nr = size(rho,1);
	ncellrho = size(area2h,1);
	N = Nr/ncellrho;


	% define functional and its derivative
	fun = fun_dfun_ddfun{1};
	der = fun_dfun_ddfun{2};
	derder = fun_dfun_ddfun{3};

	% Integration using midpoint (all orders will use middle_rho)
	rho_all = [rho_initial;rho;rho_final];
	middle_rho = cell(N+1,1);
	for k=1:N+1
    middle_rho{k} = 0.5*(...
													rho_all(1+(k-1)*ncellrho:k    *ncellrho) + ...
													rho_all(1+(k  )*ncellrho:(k+1)*ncellrho)  );
	end

	sign = 1;

	if (order == 0)
		% funcitonal
		out = 0
		for i=1:N+1
			out = out + area2h'*fun(middle_rho{k})/(N+1);
		end
	elseif ( order == 1)
		% gradient
		out = zeros(Nr,1);
		for k = 1: N
			% d F/d rho_i = 
			out(1+(k-1)*ncellrho:k*ncellrho) = 1/(N+1) * area2h .* ...
																				 (0.5*der(middle_rho{k})+...
																					0.5*der(middle_rho{k+1}));
		end
		out = sign*out; % Check this sign
	elseif (order == 2)
		% assemble the hessian
		hessian_component = cell(N+1,1);
		for k=1:N+1
			% evaluate hessian at middle rho
			diag = area2h/(N+1).*derder(middle_rho{k});
      hessian_component{k} = sparse(1:ncellrho,1:ncellrho,diag,ncellrho,ncellrho);
		end
		RHt = assembleRHt(N-1,ncellrho);
    hessian_component = blkdiag(hessian_component{:});
    out = RHt*hessian_component*RHt';
		out = sign*out; % Check this sign
	end 

end
