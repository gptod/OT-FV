function OC = Fkgeod(N,rho_f,rho_in,Dt,divt,Mxt,Mxt2h,Mst,gradt,It,rhosk,drhosk,uk,mu)

% Compute the system of first-order optimality conditions

Nt = N+1;
ncell2h = size(Mxt2h,1)/Nt;
ncell = size(Mxt,1)/Nt;
tnr2h = N*ncell2h;
tnp = Nt*ncell;
phi = uk(1:tnp);
rho = uk(tnp+1:tnp+tnr2h);
s = uk(tnp+tnr2h+1:tnp+2*tnr2h);

rho_all = [rho_in;rho;rho_f];
I_all = [sparse(ncell2h,tnr2h);speye(tnr2h,tnr2h);sparse(ncell2h,tnr2h)];
gradphi = gradt*phi;

% optimality conditions in phi
OCp = Mxt*Nt*Dt*It*rho_all + divt*(rhosk.*gradphi);

% optimality conditions in rho
OCr = (Nt*Mxt*Dt*It*I_all)'*phi - (0.5*drhosk'*Mst*(gradphi).^2) - Mxt2h(1:tnr2h,1:tnr2h)*s;

% optimality conditions in s
OCs = (-rho.*s+mu);


if (0)
  for i=1:N
    %G^k=(phi_k)^T *(F^phi_k ) + (rho_k)^T*(F^rho_k)
    pi=phi(1+(i-1)*ncell  :i*ncell  )'*OCp(1+(i-1)*ncell  :i*ncell  )+...
       rho(1+(i-1)*ncell2h:i*ncell2h)'*OCr(1+(i-1)*ncell2h:i*ncell2h);

    % 
    OCp((i+1)*ncell)=pi;    
  end

  % 1/Deltat* phi_{N+1} ^T * I * rho_N
  OCp(Np)= Nt* phi(1+N*ncell  :(N+1)*ncell  )'* ...
	      ( Mxt(1:ncell,1:ncell) * It(1:ncell,1:ncell2h)*rho(1+(N-1)*ncell2h:N*ncell2h))
  
end


OC.p = OCp;
OC.r = OCr;
OC.s = OCs;


end

