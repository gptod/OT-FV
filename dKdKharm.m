function r=dKdKharm_vect(rhoK,rhoL,dK,dL,ds)

% double derivative w.r.t. rhoK of the weighted harmonic mean

r=zeros(size(rhoK));

a = ismember([rhoK rhoL],[0 0],'rows');
b = ~a;
r(b) = -2*ds(b).*dK(b).*dL(b).*rhoL(b).^2./(dK(b).*rhoL(b)+dL(b).*rhoK(b)).^3;


end



