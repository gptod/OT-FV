function [S]=compute_SchurCA(A,B1T,B2,C,Nt,ncellrho,ncellphi)

  indc=1
  ctrl=ctrl_solver;
  ctrl.init('agmg',1e-12, 200)
  
  invAB1T=zeros(ncellrho,ncellrho,2*Nt-2)  ;
 
  Ai=A(1:ncellrho,1:ncellrho);				% ground the solution
  Ai(indc,:) = sparse(1,ncellrho);
  Ai(:,indc) = sparse(ncellrho,1);
  Ai(indc,indc) = 1;
  invAi=sparse_inverse;
  invAi.init(Ai,ctrl);

  imat=1; irow=1; icol=1;
  for j=1:ncellrho
      v=B1T( ncellrho*(irow-1)+1 : irow*ncellrho , ...
	     ncellrho*(icol-1)+j);
    % poject to ortogonal of the kernel
    rhs=v-ones(ncellrho,1)*(ones(ncellrho,1)'*v/(ncellrho));
    rhs(indc)     = 0;
    invAB1T(:,j,imat)=invAi.apply(rhs);
    invAi.info.print();
  end
  invAi.kill();

  for k=2:Nt-1
    % build A_i and its inverse
    irow=k;
    Ai=A(ncellrho*(irow-1)+1:irow*ncellrho,...
	 ncellrho*(icol-1)+1:icol*ncellrho);
    Ai(indc,:) = sparse(1,ncellrho);
    Ai(:,indc) = sparse(ncellrho,1);
    Ai(indc,indc) = 1;
    invAi=sparse_inverse;
    invAi.init(Ai,ctrl);

    % i,i-1
    icol=k-1;
    for j=1:ncellrho
      v=B1T( ncellrho*(irow-1)+1 : irow*ncellrho , ...
	     ncellrho*(icol-1)+j);
				% poject to ortogonal of the kernel
      rhs=v-ones(ncellrho,1)*(ones(ncellrho,1)'*v/(ncellrho));
      rhs(indc)     = 0;
      invAB1T(:,j,imat)=invAi.apply(rhs);
      invAi.info.print();
    end

    % i,i
    icol=k;
    for j=1:ncellrho
      v=B1T( ncellrho*(irow-1)+1 : irow*ncellrho , ...
	     ncellrho*(icol-1)+j);
				% poject to ortogonal of the kernel
      rhs=v-ones(ncellrho,1)*(ones(ncellrho,1)'*v/(ncellrho));
      rhs(indc)     = 0;
      invAB1T(:,j,imat)=invAi.apply(rhs);
      invAi.info.print();
    end

  end

				% i,i
  irow=Nt; icol=Nt-1;
  imat=imat+1
  for j=1:ncellrho
    v=B1T( ncellrho*(irow-1)+1 : irow*ncellrho , ...
	   ncellrho*(icol-1)+j);
				% poject to ortogonal of the kernel
    rhs=v-ones(ncellrho,1)*(ones(ncellrho,1)'*v/(ncellrho));
    rhs(indc)     = 0;
    invAB1T(:,j,imat)=invAi.apply(rhs);
    invAi.info.print();
  end
  


  return
