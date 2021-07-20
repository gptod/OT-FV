function [A,B1T,rhs,B1T_perturbation]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas)

  N=size(alphas,1);
  Nt=N+1;
  ncellphi=size(A,1)/Nt;
  ncellrho=size(B1T,2)/N;
  just_replace=0;

  B1T_perturbation=sparse(size(B1T,1),size(B1T,2));
  for i=2:Nt
    iblock=i-1; %1
    jblock=i;   %2
    irow=indeces_global(jblock);
    if (0)
    fprintf('%d alphas(i)=%1.4e, norm(v1)=%1.4e norm(v2)=%1.4e \n',...
	    iblock,alphas(iblock),...
	    norm(vectors_x(:,1+(iblock-1)*2)),...
	    norm(vectors_x(:,2+(iblock-1)*2)));
    end
    if (just_replace)
        A(irow,1+(jblock-2)*ncellphi:(jblock-1)*ncellphi) =sparse(1,ncellphi);
	A(irow,1+(jblock-1)*ncellphi:jblock    *ncellphi) ==sparse(1,ncellphi);
	B1T(irow,:)=sparse(1,size(B1T,2));
	rhs(irow)=0;
    end

	
    A(irow,1+(jblock-2)*ncellphi:(jblock-1)*ncellphi) = ...
    A(irow,1+(jblock-2)*ncellphi:(jblock-1)*ncellphi)+...
    vectors_x(:,1+(iblock-1)*2)';
    
    

    A(irow,1+(jblock-1)*ncellphi:jblock    *ncellphi) = ...
    A(irow,1+(jblock-1)*ncellphi:jblock    *ncellphi) + ...
    vectors_x(:,2+(iblock-1)*2)';

    
    B1T(irow,:)=B1T(irow,:)+vectors_y(:,iblock)';

    B1T_perturbation(irow,:)=vectors_y(:,iblock)';
    
    rhs(irow)=rhs(irow)+alphas(iblock);
  end


end
