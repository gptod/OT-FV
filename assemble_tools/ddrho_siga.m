function [ddrhosa] = ddrho_siga(ind,sigma,rho,a)

% second derivative of the harmonic reconstruction from cells to edges applied to
% the vector a

ncell = size(rho,1);

KK=sigma(ind.internal,3);
LL=sigma(ind.internal,4);

dK = sigma(ind.internal,5);
dL = sigma(ind.internal,6);
ds = sigma(ind.internal,3);


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




