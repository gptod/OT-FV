function [S]=assemble_dual_schur(A,B1T,B2,C,W_mat,ncellphi,ncellrho)
	Np=size(A,1);
	Nr=size(C,1);
	N=Nr/ncellrho;
	
	invA=sparse_inverse;
	ctrl_loc=ctrl_solver;
	ctrl_loc.init('agmg',...
								1e-14,...
								200,...
								0,...
								1,...
								'A',64);

	
	invA.init(A+1e-12*speye(Np,Np),ctrl_loc);
	S=zeros(Nr,Nr);
	temp=zeros(Nr,1);

	P = @(y) ort_proj(y,W_mat)
	
	for i=1:Nr
		disp(i/Nr*100);
		temp(:)=0;
		temp(i)=1;
		S(:,i)=P(C*P(temp))+P(B2*(invA.apply(B1T*P(temp))));
	end
