function [y]=prec_with_solvesys(JF,F,controls,x)
  JF.p=-x(1:Np);
  JF.p=-x(Np+1:Np+Nr);
  JF.p=-x(Np+Nr+1:Np+2*Nr);

  y=solvesys(JF,F,controls);

end
