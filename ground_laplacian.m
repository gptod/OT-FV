function A = ground_laplacian(A,inode)
  A(inode,:)=sparse(size(A,1),1)
  A(inode,inode)=1.0
end
