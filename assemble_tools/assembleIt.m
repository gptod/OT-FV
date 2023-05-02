function It = assembleIt(N,ncell,ncell2h,I)


% It = sparse((N+2)*ncell,(N+2)*ncell2h);
% 
% for k=1:N+2
%    
%     It((k-1)*ncell+1:k*ncell,(k-1)*ncell2h+1:k*ncell2h) = I;
%     
% end

It = repmat({I},1,N+2);
It = blkdiag(It{:});