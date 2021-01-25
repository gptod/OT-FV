function [resnorm] = compute_linear_system_residuum(JF,F,d)
  Np = length(F.p);
  Nr = length(F.r);
  Ns = length(F.s);
  res(1:Np)            =F.p+JF.pp*d(1:Np)      +JF.pr*d(1+Np:Np+Nr);
  res(1+Np:Np+Nr)      =F.r+JF.rp*d(1:Np)      +JF.rr*d(1+Np:Np+Nr)+JF.rs*d(Np+Nr+1:Np+Nr+Ns);
  res(Np+Nr+1:Np+Nr+Ns)=F.s+                    JF.sr*d(1+Np:Np+Nr)+JF.ss*d(Np+Nr+1:Np+Nr+Ns);

  resnorm=norm(res)/norm([F.p;F.r;F.s]);
end
