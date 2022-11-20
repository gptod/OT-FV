function [f_df_ddf] = symbolic2f_df_ddf(symbolic_function)


	rho = argnames(symbolic_function);
	
	der = diff(symbolic_function);
	derder = diff(der);

	
	f_df_ddf = {matlabFunction(symbolic_function),matlabFunction(der),matlabFunction(derder)};

	end
