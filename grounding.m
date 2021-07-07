function [A,B1T,rhs]=grounding(A,B1T,rhs,inode,value)
  Np=size(A,1);
  Nr=size(B1T,2);
  adiag=A(inode,inode);
  adiag=1.0;
  rhs(inode)=adiag*value;
  A(inode,:)   = sparse(1,Np);
  A(inode,inode)= adiag;
  B1T(inode,:) = sparse(1,Nr);

  

end
