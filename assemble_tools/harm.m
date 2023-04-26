function r=harm_vect(rhoK,rhoL,dK,dL,ds)

% weighted harmonic mean

r=zeros(size(rhoK));
a = ismember([rhoK rhoL],[0 0],'rows');
b = ~a;
r(b)= (ds(b).*rhoK(b).*rhoL(b))./(rhoL(b).*dK(b)+rhoK(b).*dL(b));