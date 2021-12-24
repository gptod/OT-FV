function [augOC] = augmented_FKgeod(OC,uk, Number_of_rows, line,integral_constrain,sign,add,verbose) 
  %				
  % Pass the ORIGINAL (before augmentation) OC 
  %
  augOC=OC;

				%
  

  N = OC.N;
  Nt=N+1;
  ncell2h=OC.ncell2h;
  ncell=OC.ncell;

  tnr2h = N*ncell2h;
  tnp = (N+1)*ncell;
  phi = uk(1:tnp);
  rho = uk(tnp+1:tnp+tnr2h);
  s = uk(tnp+tnr2h+1:tnp+2*tnr2h);

  
  for k=1:Number_of_rows
    %G^k=(phi_k)^T *(F^phi_k ) + (rho_k)^T*(F^rho_k)
    new=get_slice(phi,k,ncell,N+1)'*get_slice(OC.p,k,ncell,N+1)+...
       sign*get_slice(rho,k,ncell2h,N)'*get_slice(OC.r,k,ncell2h,N);

    % either replace or add equation 
    if ( not(add) )
      augOC.p(line+(k-1)*ncell)=0
    end
    augOC.p(line+(k-1)*ncell)=augOC.p(line+(k-1)*ncell)+new;

    % print info
    if (verbose)
      fprintf('F_phi_%d(%d) old=%1.2e new aug=%1.2e phi(line)=%1.2e |phi|=%1.2e \n', ...
	      k,line,OC.p(line+(k-1)*ncell),new,phi(k*ncell),norm(get_slice(phi,k,ncell,N+1)))
    end
  end

  if (integral_constrain)
    % G_N+1=1/Deltat* phi_{N+1} ^T * I * rho_N
    if ( not(add) )
      augOC.p(Np)=0
    end
    augOC.p(Np) = augOC.p(Np) +...
		  Nt* get_slice(phi,N+1,ncell,N+1)'* ...
		  ( Mxt(1:ncell,1:ncell) * It(1:ncell,1:ncell2h)*get_slice(rho,N,ncell2h,N));
  end
  
end
