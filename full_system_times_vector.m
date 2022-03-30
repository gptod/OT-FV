function [y] = full_system_times_vector(d,A,B1T,B2,R,M,Ds,Dr,Np,Nr)
  Ns = Nr;
  y=zeros(Np+Nr+Nr,1);
  res(1:Np,1)            =A (d(1:Np))     +B1T(d(1+Np:Np+Nr));
  res(Np+1:Np+Nr,1)      =B2(d(1:Np))     +R  (d(1+Np:Np+Nr))   +M (d(Np+Nr+1:Np+Nr+Ns));
  res(Np+Nr+1:Np+Nr+Ns,1)=                 Ds (d(1+Np:Np+Nr))   +Dr(d(Np+Nr+1:Np+Nr+Ns));
end
