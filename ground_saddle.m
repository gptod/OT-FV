function [A,B1T_time,B1T_space,B2_time,B2_space,f1]=ground_saddle(A,B1T_time,B1T_space,f1,N)
	Np=size(A,1);
  Nr=size(B1T_space,2);
	ncellphi=Np/(N+1);
	
	node=1;
	nodes=zeros(N+1,1);
	for i=1:N+1
		nodes(i)=node+(i-1)*ncellphi;
	end

	B2_time=B1T_time';
	
	B2_space=B1T_space';
		
	% columns
	if (0)
		for i=1:N+1
			indc=nodes(i);
			B2_time(:,indc)= sparse(Nr,1);
			B2_space(:,indc)= sparse(Nr,1);
			A(:,indc) = sparse(Np,1);
		end
	end
	
	% rows
	for i=1:N+1
		indc=nodes(i);
		A(indc,:) = sparse(1,Np);
		A(indc,indc)= 1e0;
		B1T_time(indc,:)  = sparse(1,Nr);
		B1T_space(indc,:) = sparse(1,Nr);
		f1(indc)=0;
	end



end
