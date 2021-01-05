function r=dLdKharm(rhoK,rhoL,dK,dL,ds)

r = 2*ds.*dK.*dL.*rhoK.*rhoL./(dK.*rhoL+dL.*rhoK).^3;

end