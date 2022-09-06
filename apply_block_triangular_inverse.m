function [sol]=apply_block_triangular_inverse(diagonal_inverse,extra_diagonal_matrices,L_or_U,rhs)
  nblocks=size(diagonal_inverse,1);
  nrow=size(rhs,1);
  nrow_block=size(diagonal_inverse(nblocks).matrix,1);
  %fprintf('%d %d \n', nblocks, nrow)
  sol=zeros(nrow,1);
  %for i=1,nblocks
  %  diagonal_inverse(i).info()
  %end
  if ( strcmp(L_or_U,'L') ) 
    temp=rhs;    
    nrow_block=size(diagonal_inverse(1).matrix,1);
    ibegin = 1;
    iend   = nrow_block;
    ibegin_before = ibegin;
    iend_before   = iend;
    sol(ibegin:iend) = diagonal_inverse(1).apply(temp(ibegin:iend));
    %fprintf('%d %d \n', ibegin, iend)
    for i=2:nblocks
      nrow_block=size(diagonal_inverse(i).matrix,1);
      ibegin = ibegin_before+ nrow_block;
      iend   = iend_before  + nrow_block;
      %fprintf('%d %d | %d %d \n', ibegin, iend,ibegin_before,iend_before)
      temp(ibegin:iend) = temp(ibegin:iend)...
			  - extra_diagonal_matrices{(i-1)}*sol(ibegin_before:iend_before);
      sol(ibegin:iend)=diagonal_inverse(i).apply(temp(ibegin:iend));
      ibegin_before = ibegin;
      iend_before   = iend;
    end
    
  elseif  ( strcmp(L_or_U,'U') )
    temp=rhs;
    nrow_block=size(diagonal_inverse(nblocks).matrix,1);

    
    ibegin= nrow-nrow_block+1;
    iend  = nrow;
    ibegin_before = ibegin;
    iend_before   = iend;
    
    sol(ibegin:iend) = diagonal_inverse(nblocks).apply(temp(ibegin:iend));
    
    for j=1:nblocks-1
      i=nblocks-j;
      nrow_block=size(diagonal_inverse(i).matrix,1);
      ibegin = ibegin_before-nrow_block;
      iend   = iend_before-nrow_block;
      %fprintf('%d %d %d %d %d  \n',i,ibegin,iend,ibegin_before,iend_before) 
      temp(ibegin:iend) = temp(ibegin:iend) - ...
			  extra_diagonal_matrices{i}*sol(ibegin_before:iend_before);
      sol(ibegin:iend)=diagonal_inverse(i).apply(temp(ibegin:iend));
      ibegin_before = ibegin;
      iend_before   = iend;
    end
  else
    disp('String can be L or U')
    return
  end
