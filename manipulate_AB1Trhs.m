function [A,B1T,rhs]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas)

  N=size(alphas,1);
  Nt=N+1;
  ncellphi=size(A,1)/Nt;
  ncellrho=size(B1T,2)/N;
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
    %A(irow,1+(jblock-2)*ncellphi:(jblock-1)*ncellphi) =sparse(1,ncellphi);
    A(irow,1+(jblock-2)*ncellphi:(jblock-1)*ncellphi) = ...
    A(irow,1+(jblock-2)*ncellphi:(jblock-1)*ncellphi)+...
    vectors_x(:,1+(iblock-1)*2)';
    
    %A(irow,1+(jblock-1)*ncellphi:jblock    *ncellphi) ==sparse(1,ncellphi);
    A(irow,1+(jblock-1)*ncellphi:jblock    *ncellphi) = ...
    A(irow,1+(jblock-1)*ncellphi:jblock    *ncellphi) + ...
    vectors_x(:,2+(iblock-1)*2)';

    %B1T(irow,:)=sparse(1,size(B1T,2));
    B1T(irow,:)=B1T(irow,:)+vectors_y(:,iblock)';

    %rhs(irow)=0;
    rhs(irow)=rhs(irow)+alphas(iblock);
  end
end
