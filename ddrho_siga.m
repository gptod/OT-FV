function [ddrhosa] = ddrho_siga(ind,edges,mid,cc,rho,a)

% second derivative of the harmonic reconstruction from cells to edges applied to
% the vector a

ncell = size(cc,1);
nei = length(ind.internal);

% for e=1:nei
%     
%     K = edges(ind.internal(e),3);
%     L = edges(ind.internal(e),4);
%     dK = sqrt((mid(ind.internal(e),1)-cc(K,1))^2+(mid(ind.internal(e),2)-cc(K,2))^2);
%     dL = sqrt((mid(ind.internal(e),1)-cc(L,1))^2+(mid(ind.internal(e),2)-cc(L,2))^2);
%     ds = edges(ind.internal(e),5);
%     % harmonic mean
%     ddrhosa(K,K) = ddrhosa(K,K) + dKdKharm(rho(K),rho(L),dK,dL,ds)*a(e);
%     ddrhosa(K,L) = ddrhosa(K,L) + dLdKharm(rho(K),rho(L),dK,dL,ds)*a(e);
%     ddrhosa(L,K) = ddrhosa(L,K) + dLdKharm(rho(K),rho(L),dK,dL,ds)*a(e);
%     ddrhosa(L,L) = ddrhosa(L,L) + dKdKharm(rho(L),rho(K),dL,dK,ds)*a(e);
% 
% end


KK=edges(ind.internal,3);
LL=edges(ind.internal,4);

% % possono essere generati nella struttura ind
% indKK = sub2ind([ncell ncell],KK,KK);
% indKL = sub2ind([ncell ncell],KK,LL);
% indLK = sub2ind([ncell ncell],LL,KK);
% indLL = sub2ind([ncell ncell],LL,LL);
% [indu_KK,~,indo_KK] = unique(indKK);
% [indu_KL,~,indo_KL] = unique(indKL);
% [indu_LK,~,indo_LK] = unique(indLK);
% [indu_LL,~,indo_LL] = unique(indLL);


% in futuro aggiungiamo il calcolo di dK e dL in edges
% dK = sqrt((mid(ind.internal(:),1)-cc(KK,1)).^2+...
%     (mid(ind.internal(:),2)-cc(KK,2)).^2);
% dL = sqrt((mid(ind.internal(:),1)-cc(LL,1)).^2+...
%     (mid(ind.internal(:),2)-cc(LL,2)).^2);
dK = edges(ind.internal,8);
dL = edges(ind.internal,9);
ds = edges(ind.internal,5);


% ddrhosa = sparse(ncell,ncell);
% ddrhosa(ind.u_KK) = ddrhosa(ind.u_KK) + ...
%                    accumarray(ind.o_KK,dKdKharm(rho(KK),rho(LL),dK,dL,ds).*a);
% ddrhosa(ind.u_KL) = ddrhosa(ind.u_KL) + ...
%                    accumarray(ind.o_KL,dLdKharm(rho(KK),rho(LL),dK,dL,ds).*a);         
% ddrhosa(ind.u_LK) = ddrhosa(ind.u_LK) + ...
%                    accumarray(ind.o_LK,dLdKharm(rho(KK),rho(LL),dK,dL,ds).*a);          
% ddrhosa(ind.u_LL) = ddrhosa(ind.u_LL) + ...
%                    accumarray(ind.o_LL,dKdKharm(rho(LL),rho(KK),dL,dK,ds).*a);           


% no need to use the accumarray function as the sparse function sums by
% default values with the same indices
% this assigning is not faster than the other one
ddrhosa = sparse(KK,KK,dKdKharm(rho(KK),rho(LL),dK,dL,ds).*a,ncell,ncell);
ddrhosa = ddrhosa + sparse(KK,LL,dLdKharm(rho(KK),rho(LL),dK,dL,ds).*a,ncell,ncell);
ddrhosa = ddrhosa + sparse(LL,KK,dLdKharm(rho(KK),rho(LL),dK,dL,ds).*a,ncell,ncell);
ddrhosa = ddrhosa + sparse(LL,LL,dKdKharm(rho(LL),rho(KK),dL,dK,ds).*a,ncell,ncell);




