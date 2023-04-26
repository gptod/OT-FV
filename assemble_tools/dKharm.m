function r=dKharm_vect(rhoK,rhoL,dK,dL,ds)

% derivative w.r.t. rhoK of the weighted harmonic mean

r=zeros(size(rhoK));
a = rhoL~=0;
r(a)= ds(a).*dK(a)./(dK(a)+dL(a).*rhoK(a)./rhoL(a)).^2;

end