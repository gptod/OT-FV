function [sol]=apply_inverseAtilde(array_sparse_inverse,vectors,indeces_local,rhs)
  nblocks=size(array_sparse_inverse,1);
  start=1;
  sol=zeros(size(rhs,1),1);
  %disp('Block inverse A')

  %start
  nrow_block=size(array_sparse_inverse(1).matrix,1);
  sol(1:nrow_block)=array_sparse_inverse(1).apply(rhs(1:nrow_block));
  start=nrow_block+1;
  for i=2:nblocks    
    %start
    nrow_block=size(array_sparse_inverse(i).matrix,1);
    vec=rhs(start:start+nrow_block-1);
    vec(indeces_local(i))=vec(indeces_local(i))-vectors(:,1+(i-2)*2)'*sol(start-nrow_block:start-1);
    sol(start:start+nrow_block-1)=array_sparse_inverse(i).apply(vec);
    start=start+nrow_block;
  end
