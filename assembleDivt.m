function Divt = assembleDivt(N,ncell,nedge,div)


% Divt = sparse((N+1)*ncell,(N+1)*nedge);
% 
% for k=1:N+1
%    
%     Divt((k-1)*ncell+1:k*ncell,(k-1)*nedge+1:k*nedge) = div;
%     
% end

Divt = repmat({div},1,N+1);
Divt = blkdiag(Divt{:});