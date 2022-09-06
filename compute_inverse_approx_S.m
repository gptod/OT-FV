function [x]=compute_inverse_approx_S(y, inverse,A_rho,pr,lrb,mode_Arho,invC,Dr,Mrho,inv_Mrho)
	if (mode_Arho==1)
		if strcmp(lrb,'l')
			x=inverse(A_rho*pr'*pr*y);
		elseif strcmp(lrb,'r')
			x=pr'*pr*A_rho*(inverse(y));
		elseif strcmp(lrb,'b')
			x=pr'*pr*A_rho*(inverse(A_rho*pr'*pr*y));	
		end
	elseif (mode_Arho==2)
		if strcmp(lrb,'l')
			x=inverse(A_rho*y);
		elseif strcmp(lrb,'r')
			x=inv_Mrho*A_rho*(inverse(Dr*inv_Mrho*y));
		elseif strcmp(lrb,'b')
			x=A_rho*(inverse(A_rho*y));
		end
	end
		
