function r=dLdKharm_vect(rhoK,rhoL,dK,dL,ds)

% double derivative w.r.t. rhoK, rhoL of the weighted harmonic mean

r=zeros(size(rhoK));

a = ismember([rhoK rhoL],[0 0],'rows');
b = ~a;
r(b) = 2*ds(b).*dK(b).*dL(b).*rhoK(b).*rhoL(b)./(dK(b).*rhoL(b)+dL(b).*rhoK(b)).^3;


end