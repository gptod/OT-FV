function r=dKharm(rhoK,rhoL,dK,dL,ds)

r = ds.*dK./(dK+dL.*rhoK./rhoL).^2;

end