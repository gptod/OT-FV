function [res] = J_times_vectors(JF,d)
  Np = size(JF.pp,1);
  Nr = size(JF.rr,1);
  Ns = size(JF.ss,1);
  res=zeros(Np+Nr+Ns,1);
  res(1:Np)            =JF.pp*d(1:Np)      +JF.pr*d(1+Np:Np+Nr);
  res(1+Np:Np+Nr)      =JF.rp*d(1:Np)      +JF.rr*d(1+Np:Np+Nr)+JF.rs*d(Np+Nr+1:Np+Nr+Ns);
  res(Np+Nr+1:Np+Nr+Ns)=                    JF.sr*d(1+Np:Np+Nr)+JF.ss*d(Np+Nr+1:Np+Nr+Ns);
end
