function [A,B1T,B2,rhs]=grounding(A,B1T,B2,rhs,inode,value)
  Np=size(A,1);
  Nr=size(B1T,2);

  %A(:,inode)=sparse(Np,1);
  %B2(:,inode)=sparse(Nr,1);
  

  

  B1T(inode,:) = sparse(1,Nr);

  A(inode,:)   = sparse(1,Np);
  A(inode,inode)= 1.0;
  rhs(inode)=value;
end
