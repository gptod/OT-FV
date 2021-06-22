function gradt = assembleGradt(N,ncell,nedge,grad)

% gradt = sparse((N+1)*nedge,(N+1)*ncell);
% 
% for k=1:N+1
%    
%     gradt((k-1)*nedge+1:k*nedge,(k-1)*ncell+1:k*ncell) = grad;
%     
% end

gradt = repmat({grad},1,N+1);
gradt = blkdiag(gradt{:});