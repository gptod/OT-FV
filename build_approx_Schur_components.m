function [pr,A_rho, approx_S,start]=build_approx_Schur_components(A,B1T,B2,C,lrb,JF,mode,mode_Arho)
	% lrb = l=left,r=right,b=both 

	
	% get dimensions
	ncellrho = JF.ncellrho;	
	ncellphi = JF.ncellphi;
	N=JF.ntimestep;
	Np=size(A,1);
	Nr=size(B2,1);
	
	
	% define operator pr that project rho_variable into phi_variables
	I=JF.I;
	Pspace=assemble_space_projector(N+1,I');  
	Ptime=assemble_time_projector(N,ncellrho);

	% weights
	Mphi=JF.Mxt;
	Mrho=-JF.rs;
	inv_Mphi  = sparse(1:Np,1:Np,(1./spdiags(Mphi,0))',Np,Np);
	inv_Mrho  = sparse(1:Nr,1:Nr,(1./spdiags(Mrho,0))',Nr,Nr);
		
	
	
	pr = Pspace'*Ptime';

	if (mode_Arho==1)
		A_rho = Ptime*Pspace*inv_Mphi*A*inv_Mphi*Pspace'*Ptime';
	elseif (mode_Arho==2)		
		old=0;
		if (old)
			nei=size(JF.div,2);
			divt = assembleDivt(N-1,ncellrho,nei,JF.div);
			gradt = assembleGradt(N-1,ncellrho,nei,JF.grad);
			neit=nei*N;
			Rst = repmat({JF.Rs},1,N);
			Rst = blkdiag(Rst{:});
			rho=spdiags(-JF.ss,0);
			rho_edge=spdiags(Rst*rho,0,neit,neit);
			if (strcmp(mode,'commute'))
				A_rho = -inv_Mrho * divt*(rho_edge)*gradt*inv_Mrho;
			elseif (strcmp(mode,'commute_bis'))
				A_rho = -divt*(rho_edge)*gradt*inv_Mrho;
			elseif (strcmp(mode,'commute_ter'))
				A_rho = - divt*(rho_edge)*gradt;
			end
		else		
			neit=size(JF.divt_rho,2);
			rho=spdiags(-JF.ss,0);
			rho_edge=spdiags(JF.Rst_rho*rho,0,neit,neit);
			%A_rho = -inv_Mrho * JF.divt_rho*(rho_edge)*JF.gradt_rho*inv_Mrho ;
			if (strcmp(mode,'commute'))
				A_rho = -inv_Mrho *  JF.divt_rho*(rho_edge)*JF.gradt_rho*inv_Mrho;
			elseif (strcmp(mode,'commute_bis'))
				A_rho = - JF.divt_rho*(rho_edge)*JF.gradt_rho*inv_Mrho;
			elseif (strcmp(mode,'commute_ter'))
				A_rho = -  JF.divt_rho*(rho_edge)*JF.gradt_rho;
			elseif (strcmp(mode,'commute_quater'))
				A_rho = - inv_Mrho*JF.divt_rho*(rho_edge)*JF.gradt_rho;
			end
		end
		%size(JF.divt_rho)
		%size(rho_edge)
		%size(JF.gradt_rho)
		%size(inv_Mrho)
	
	end

	start=0;

	%approx_S(start+1:start+Nr,start+1:start+Nr) = JF.B1T_time'*inv_Mphi*inv_Mphi*JF.B1T_time;
	%approx_S(start+1:start+Nr,start+1:start+Nr) = B2*inv_Mphi*B1T;
	%approx_S(start+1:start+Nr,start+1:start+Nr) = JF.B1T_time'*inv_Mphi*inv_Mphi*JF.B1T_time+...
	%																							JF.B1T_space'*inv_Mphi*inv_Mphi*JF.B1T_time;
	
	if strcmp(lrb,'l')
		approx_S(start+1:start+Nr,start+1:start+Nr) = B2*inv_Mphi*inv_Mphi*B1T;	
		approx_S(start+1:start+Nr,start+1:start+Nr) = A_rho*pr'*pr*C+approx_S(start+1:start+Nr,start+1:start+Nr);
	elseif strcmp(lrb,'r')
		if (strcmp(mode,'commute'))
			approx_S(start+1:start+Nr,start+1:start+Nr) = B2*inv_Mphi*inv_Mphi*B1T;	
		elseif (strcmp(mode,'commute_bis'))
			approx_S(start+1:start+Nr,start+1:start+Nr) = B2*inv_Mphi*B1T;	
		elseif (strcmp(mode,'commute_ter'))
			approx_S(start+1:start+Nr,start+1:start+Nr) = B2*B1T;
		elseif (strcmp(mode,'commute_quater'))
			approx_S(start+1:start+Nr,start+1:start+Nr) = B2*inv_Mphi*B1T;	
		end		
		%approx_S(start+1:start+Nr,start+1:start+Nr) = C*A_rho*pr'*pr+approx_S(start+1:start+Nr,start+1:start+Nr);
		approx_S(start+1:start+Nr,start+1:start+Nr) = C*A_rho+approx_S(start+1:start+Nr,start+1:start+Nr);
	elseif strcmp(lrb,'b')
		disp('here')
		%approx_S(start+1:start+Nr,start+1:start+Nr) = JF.B1T_time'*inv_Mphi*inv_Mphi*A*inv_Mphi*inv_Mphi*JF.B1T_time;
		approx_S(start+1:start+Nr,start+1:start+Nr) = B2*inv_Mphi*inv_Mphi*A*inv_Mphi*inv_Mphi*B1T;	
		approx_S(start+1:start+Nr,start+1:start+Nr) = approx_S(start+1:start+Nr,start+1:start+Nr) + ...
																									pr'*pr*A_rho*C*A_rho*pr'*pr;
																									
	end

	% %return
	% if strcmp(lrb,'l')
	% 	approx_S = A_rho*pr'*pr*C+JF.Dt'*inv_Mphi*inv_Mphi*JF.Dt;
	% elseif strcmp(lrb,'r')
	% 	approx_S = C*A_rho*pr'*pr+JF.Dt'*inv_Mphi*inv_Mphi*JF.Dt;
	% elseif strcmp(lrb,'b')
	% 	approx_S = pr'*pr*A_rho*C*A_rho*pr'*pr+JF.Dt'*inv_Mphi*inv_Mphi*A*inv_Mphi*inv_Mphi*JF.Dt;
	% end


