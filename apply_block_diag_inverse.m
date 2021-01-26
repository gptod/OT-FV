%%

%% matrices = cell(4, 1);
%% for i = 1:4    
%%  A{i} = speye(4);
%% end
%% to init and array of inverse
%% matrices(1:nblock)=speye(4)
%% array_inv(1:nblock)=sparse_inverse
%% for i=1:nblock
%%    array_inv(i).init(matrices(i),ctrl)
%% end 
%%

function [x]=apply_block_diag_inverse(array_sparse_inverse,rhs)
  nblocks=size(array_sparse_inverse,1)
  start=1
  for i=1:nblocks
    nrow_block=size(array_sparse_inverse(i).matrix,1)
    y(start:nrow_block)=array_sparse_inverse(i).apply(rhs(start:nrow_block))
    start=start+nrow_block
  end
