function [A,B1T,B2,C,f1,f2,W_mat]=get_saddle_point(JF,F)
	N = JF.ntimestep;
	Nr = size(JF.rr,1);
	ncellrho = JF.ncellrho;
	

	% rhs
	f=-F.p;
	g=-F.r;
	h=-F.s;

	
	% Solve with the Schur complement A|C, using Matlab's backslash (or agmg)
	swap_sign=1;
	sign=-1;
	
	% reduced system to phi rho variables
	A = sparse(JF.pp);
	B1T = sparse(JF.pr-JF.ps*(JF.ss\JF.sr));
	B2 = sparse(JF.rp);
	R = sparse(JF.rr);
	M   = sparse(JF.rs);
	Ds = sparse(JF.sr);
	Dr = sparse(JF.ss);
	C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
	temp=(JF.ss\h);
	f1 = f-JF.ps*temp;
	f2 = g-JF.rs*temp;

	% swap C sign for having standard saddle point notation 
	C = -C;
	if (swap_sign)
		A=-A;
		B1T=-B1T;
		B2=-B2;
		C=-C;
		
		R=-R;
		M=-M;
		Ds=-Ds;
		Dr=-Dr;

		
		f=-f;
		g=-g;
		h=-h;
		
		f1=-f1;
		f2=-f2;
		
	end

	W_mat=zeros(N,Nr);
	for k = 1:N
	  W_mat(k,1+(k-1)*ncellrho:k*ncellrho) = JF.area2h';
  end
