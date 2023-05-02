function [sqrt_diagA,sqrt_diagC,inv_sqrt_diagA,inv_sqrt_diagC] = get_diagonal_scaling (A,C)
	Np=size(A,1);
  Nr=size(C,1);


	if (0)
    fprintf('(%9.4e <= A <= %18.12e) \n',min(spdiags(A,0)),max(spdiags(A,0)));
    fprintf('(%9.4e <= C <= %18.12e) \n',min(spdiags(C)),max(spdiags(C)));
  end

	sqrt_A=sqrt(spdiags(A,0));
	sqrt_C=sqrt(spdiags(C,0));

	inv_sqrt_diagA = sparse(1:Np,1:Np,(1.0./sqrt_A)',Np,Np);
  inv_sqrt_diagC = sparse(1:Nr,1:Nr,(1.0./sqrt_C)',Nr,Nr);

	
	sqrt_diagA=sparse(1:Np,1:Np,sqrt_A',Np,Np);
	sqrt_diagC=sparse(1:Nr,1:Nr,sqrt_C',Nr,Nr);
	
	
