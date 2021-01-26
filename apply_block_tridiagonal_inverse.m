function [x]=apply_block_tridiagonal_inverse(diagonal_inverse,extra_diagonal_matrices,L_or_U,rhs)
  nblocks=size(array_sparse_inverse,1)
  if ( strcmp(L_or_U,'L') ) 
    nrow=sum(size(array_sparse_inverse(:).matrix,1))
    temp=rhs
    
    start=1
    nrow_block=size(array_sparse_inverse(1).matrix,1)
    ibegin= start
    iend  = nrow_block
    ibegin_before = ibegin
    iend_before   = iend
    y(ibegin:iend) = array_sparse_inverse(1).apply(temp(ibegin:iend))
    
    for i=2:nblocks
      nrow_block=size(array_sparse_inverse(i).matrix,1)
      ibegin = ibegin_before-nrow_block
      iend   = iend_before-nrow_block
      temp(ibegin:iend) = temp(ibegin:iend) - ...
			  extra_diagonal_matrices{i}*y(ibegin_before:iend_before)
      y(start:start+nrow_block)=array_sparse_inverse(i).apply(temp(start:start+nrow_block))
      ibegin_before = ibegin
      iend_before   = iend
    end
    
  elseif  ( strcmp(L_or_U,'U') )
    nrow=sum(size(array_sparse_inverse(:).matrix,1))
    temp=rhs
    
    start=nrow-size(array_sparse_inverse(nblocks).matrix,1)
    nrow_block=size(array_sparse_inverse(nblocks).matrix,1)
    ibegin= start
    iend  = nrow
    ibegin_before = ibegin
    iend_before   = iend
    y(ibegin:iend) = array_sparse_inverse(nblocks).apply(temp(ibegin:iend))
    
    for i=nblocks-1:1
      nrow_block=size(array_sparse_inverse(i).matrix,1)
      ibegin = ibegin_before-nrow_block
      iend   = iend_before-nrow_block
      temp(ibegin:iend) = temp(ibegin:iend) - ...
			  extra_diagonal_matrices{i}*y(ibegin_before:iend_before)
      y(start:start+nrow_block)=array_sparse_inverse(i).apply(temp(start:start+nrow_block))
      ibegin_before = ibegin
      iend_before   = iend
    end
  else
    disp('String can be L or U')
    return
  end
