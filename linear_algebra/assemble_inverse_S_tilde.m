function [inverse_S] = assemble_inverse_S_tilde(N,I,A,B1T_time,B1T_space,C,Mxt,P)
	Nr=size(C,1);
	Np=size(A,1);
	Nt=N+1;
	ncellrho=Nr/N;
	ncellphi=Np/Nt;

	% 
	Dt = assembleDt(N,ncellphi);
	I_all = [sparse(ncellrho,Nr);speye(Nr,Nr);sparse(ncellrho,Nr)];
	It = assembleIt(N,ncellphi,ncellrho,I);
	
	% define operator that project rho_variable into phi_variables
	Pspace=assemble_space_projector(N+1,I');  
	Ptime=assemble_time_projector(N,ncellrho);
	P_pr = Pspace'*Ptime';



	% built time laplacian
	time_Laplacian=(Nt*Dt*It*I_all)'*Mxt*(Nt*Dt*It*I_all);
	
	A_rho=P_pr'*A*Mxt*P_pr;
	approx=C*A_rho + time_Laplacian;

	% define the inverse of approx=(C*L_rho+L_time)
	inv_approx=sparse_inverse;
	ctrl_loc=ctrl_solver;
	ctrl_loc.init('agmg',...
								1e-14,...
								200,...
								0,...
								0,...
								'approx',65);
	inv_approx.init(approx+1e-12*speye(Nr,Nr),ctrl_loc);

	%
	inverse_S = @(y) P(...% project to zero avg
											P_pr'*... % project to rho space
											inv_A.apply(...%invert laplacians
																	 Mxt*...% integrate rhs
																	 P_pr*... % interpolate of to phi space
																	 inv_approx.apply(...% apply (C*L_rho+L_time)^{-1}
																										 P(y)...%project to zero avg
																									 )));
	

end
