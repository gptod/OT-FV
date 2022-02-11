function [y]=apply_Schur_inverse(x,inverse_S,P,S)


	%fprintf('imbalance x=%1.1e \n ', norm(x-P(x)))
	y=inverse_S(x);
	%fprintf('error schur=%1.1e \n ', norm(P(S*P(y))-x)/norm(x))
									
