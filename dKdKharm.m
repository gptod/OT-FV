function r=dKdKharm(rhoK,rhoL,dK,dL,ds)

r = -2*ds.*dK.*dL.*rhoL.^2./(dK.*rhoL+dL.*rhoK).^3;

end