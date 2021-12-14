function [augOC] = augmented_FKgeod(OC,uk,Nof) 
  %				
  % Pass the ORIGINAL (before augmentation) OC 
  %
  augOC=OC;

  N = OC.N;
  Nt=N+1;
  ncell2h=OC.ncell2h;
  ncell=OC.ncell;

  tnr2h = N*ncell2h;
  tnp = (N+1)*ncell;
  phi = uk(1:tnp);
  rho = uk(tnp+1:tnp+tnr2h);
  s = uk(tnp+tnr2h+1:tnp+2*tnr2h);

  if (nargin == 2)
    Number_of_rows=N;
  else
    Number_of_rows=Nof;
  end
  
  for k=1:Number_of_rows
    %G^k=(phi_k)^T *(F^phi_k ) + (rho_k)^T*(F^rho_k)
    new=get_slice(phi,k,ncell,N+1)'*get_slice(OC.p,k,ncell,N+1)+...
       get_slice(rho,k,ncell2h,N)'*get_slice(OC.r,k,ncell2h,N);
    
 
    % p1=phi1^T Fp1 + rho1^T Fr1 -> last row first block
    % p1=phi2^T Fp1 + rho2^T Fr2 -> last row second block 
    augOC.p(k*ncell)=new;
    fprintf('%d F_phi=%f phi=%f aug=%f\n', k,OC.p(k*ncell),phi(k*ncell),new)
  end

  if (0)
    % G_N+1=1/Deltat* phi_{N+1} ^T * I * rho_N
    augOC.p(Np)= Nt* get_slice(phi,N+1,ncell,N+1)'* ...
		 ( Mxt(1:ncell,1:ncell) * It(1:ncell,1:ncell2h)*get_slice(rho,N,ncell2h,N));
  end
  
end
