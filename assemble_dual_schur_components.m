function [S_tt,S_tx,S_xx]=assemble_dual_schur_components(A,Nr,B1T_time,B1T_space,B1_time,B1_space)
	Np=size(A,1);
	invA=sparse_inverse;
	ctrl_loc=ctrl_solver;
	ctrl_loc.init('agmg',...
								1e-14,...
								200,...
								0,...
								0,...
								'A',64);

	invA.init(A+1e-12*speye(Np,Np),ctrl_loc);
	S_tt=zeros(Nr,Nr);
	S_tx=zeros(Nr,Nr);
	S_xx=zeros(Nr,Nr);

	temp=zeros(Nr,1);

	
	for i=1:Nr
		disp(i/Nr*100);
		temp(:)=0;
		temp(i)=1;
		vt=invA.apply(B1T_time(temp));
		vx=invA.apply(B1T_space(temp));
		
		S_tt(:,i)=B1_time(vt);
		S_tx(:,i)=B1_time(vx);
		S_xx(:,i)=B1_space(vx);
	end
	invA.kill()
