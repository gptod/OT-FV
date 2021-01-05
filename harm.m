function r=harm(rhoK,rhoL,dK,dL,ds)

% weighted harmonic mean

if rhoK==0 && rhoL==0
    r = 0;
else   
    r = (ds.*rhoK.*rhoL)./(rhoL.*dK+rhoK.*dL);
end
