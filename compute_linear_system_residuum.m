function [resnorm,resp,resr,ress] = compute_linear_system_residuum(JF,F,d)
  Np = length(F.p);
  Nr = length(F.r);
  Ns = length(F.s);

  rhs=-[F.p;F.r;F.s];
  res=J_times_vectors(JF,d)-rhs;
  resnorm=norm(res)/norm(rhs);

  
  resp=norm(res(1:Np))/norm(F.p);
  resr=norm(res(1+Np:Np+Nr))/norm(F.r);
  ress=norm(res(Np+Nr+1:Np+Nr+Ns))/norm(F.s);
end
