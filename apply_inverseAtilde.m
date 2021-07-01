function [sol]=apply_inverseAtilde(array_sparse_inverse,Amatrix,vectors_x,indeces_global,indeces_local,rhs)
  nblocks=size(array_sparse_inverse,1);
  start=1;
  sol=zeros(size(rhs,1),1);
				%disp('Block inverse A')
  size(rhs);

  %start
  nrow_block=size(array_sparse_inverse(1).matrix,1);
  ncol_block=size(array_sparse_inverse(1).matrix,2);
  
  total_ncol=ncol_block;
  
  sol(1:nrow_block)=array_sparse_inverse(1).apply(rhs(1:nrow_block));
  start=nrow_block+1;
  for i=2:nblocks    
    %start
    nrow_block=size(array_sparse_inverse(i).matrix,1);
    ncol_block=size(array_sparse_inverse(i).matrix,2);

    
    total_ncol=total_ncol+ncol_block;
    total_ncol-ncol_block+1;
    
    % local rhs
    local_rhs=rhs(start:start+nrow_block-1);
    % taking the vector p1 taking the indeces_global(i)-row of block i,i-1
    vec1=full(Amatrix(indeces_global(i),total_ncol-2*ncol_block+1:total_ncol-ncol_block))';
    %norm(vec1);
    %j=1+2*(i-2);
    %diff=norm(vec1-vectors_x(:,j));
    %vec1=vectors_x(:,j);
    

    % change rhs at local row with p1_i^x_i
    local_rhs(indeces_local(i))=local_rhs(indeces_local(i))-vec1'*sol(start-nrow_block:start-1);

    % solsve A_i x_i = local_rhs
    sol(start:start+nrow_block-1)=array_sparse_inverse(i).apply(local_rhs);
    start=start+nrow_block;
  end

  % exact=Amatrix\rhs;

  % res=norm(Amatrix*exact-rhs)/norm(rhs);
  % res=norm(Amatrix*sol-rhs)/norm(rhs);

  % norm(sol-exact)/norm(rhs)
