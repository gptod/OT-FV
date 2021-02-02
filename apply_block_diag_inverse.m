function [sol]=apply_block_diag_inverse(array_sparse_inverse,rhs)
  nblocks=size(array_sparse_inverse,1);
  start=1;
  sol=zeros(size(rhs,1),1);
  for i=1:nblocks
    %start
    nrow_block=size(array_sparse_inverse(i).matrix,1);
    %start+nrow_block-1
    sol(start:start+nrow_block-1)=array_sparse_inverse(i).apply(rhs(start:start+nrow_block-1));
    start=start+nrow_block;
  end
