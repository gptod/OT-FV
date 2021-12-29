function [OC,JOC] = newton_augmentation(uk,OC,JOC,option,line)

  N=JOC.ntimestep;
  ncell=JOC.ncellphi;
  ncell2h=JOC.ncellrho;

  add=1;
  
  nu=size(uk,1);
  Nr = N*ncell2h;
  Np = (N+1)*ncell;
  phi = uk(1      :Np     );
  rho = uk(Np+1   :Np+Nr  );
  s   = uk(Np+Nr+1:Np+2*Nr);

  % define the aug_matrix M (called aug_matrix) multiplying the non-linear system
  % G=M * OC
  % Define also the aug_matrix P (called tensor_times_OC) such that
  % 
  % JG = M * JOC  + tensor(dM/du)* OC = M * JOC + P  
  
  if (option==1)
    % M=I + P
    % P=(0 N 0)
    %   (0 0 0)
    %   (0 0 0)
    % N^k = Id +(0             ) 
    %           (  0           )
    %  line     (1 1 1 .... 1 1)   
    %           (           0  )
    %           (             0)
    %  
    %
    %
    aug_matrix=speye(nu);
    for k=1:N
      row=line+(k-1)*ncell;
      lcol=Np+(k-1)*ncell2h+1;
      rcol=Np+(k  )*ncell2h;
      if (not(add))
	aug_matrix(row,lcol:rcol)=sparse(1:ncell2h);
      end
      aug_matrix(row,lcol:rcol)=aug_matrix(line,lcol:rcol)+1e-1*ones(1,ncell2h);
    end
    tensor_times_OC=sparse(nu,nu);
  end


  if (option==2)
    deltat=1/N;
    % M^k = ( 1                  )         ( 1                  )       
    %       (  1                 )         (  1                 )
    %  line (phi_k^1 .... phi_k^N)   ..... (rho_k^1 .... rho_k^N 
    %       (                    )         (                 1  )
    %       (                   1)         (                   1)
    %  
    %
    %
    aug_matrix=speye(nu);
    tensor_times_OC=sparse(nu,nu);

    for k=1:N
      row =line+(k-1)*ncell;
      lcol=(k-1)*ncell+1;
      rcol=k*ncell;
      if (not(add))
	aug_matrix(row,lcol:rcol)=sparse(1:ncell);
      end
      aug_matrix(row,lcol:rcol) =aug_matrix(row,lcol:rcol)+...
				 deltat*get_slice(phi,k,ncell,N+1)'; %v()=phi_k  vector depends on phi
      % d vector/ d phi_k * Fphi_k
      tensor_times_OC(row,lcol:rcol)=deltat*speye(ncell)*get_slice(OC.p,k,ncell,N+1);
    end

    for k=1:N
      row=line+(k-1)*ncell;
      lcol=Np+(k-1)*ncell2h+1;
      rcol=Np+k*ncell2h;

      if (not(add))
	aug_matrix(row,lcol:rcol)=sparse(1:ncell2h);
      end
      aug_matrix(row,lcol:rcol) = aug_matrix(row,lcol:rcol) + ...
				  deltat*(-1.0)*get_slice(rho,k,ncell2h,N)'; %v()=rho_k  vector depends on phi
      % d vector/ d rho_k * Frho_k
      tensor_times_OC(row,lcol:rcol)=-deltat*speye(ncell2h)*get_slice(OC.r,k,ncell2h,N);
    end
    
      end
    
  % augemented non-linear funciton
  OCvector=[OC.p;OC.r;OC.s];
  G=aug_matrix*OCvector;
  OC.p=G(1:Np);
  OC.r=G(1+Np:Np+Nr);
  OC.s=G(1+Np+Nr:Np+Nr+Nr);

  % augemted jacobian
  original=sparse(...
		   [JOC.pp, JOC.pr, JOC.ps;
		    JOC.rp, JOC.rr, JOC.rs;
		    JOC.rp, JOC.sr, JOC.ss]...
		 );
  JOCaug_matrix=aug_matrix*original+tensor_times_OC;
  
  JOC.pp=sparse(JOCaug_matrix(1:Np,1      :Np      ));  
  JOC.pr=sparse(JOCaug_matrix(1:Np,1+Np   :Np+Nr   ));
  JOC.ps=sparse(JOCaug_matrix(1:Np,1+Np+Nr:Np+Nr+Nr));
  
  
  JOC.rp=sparse(JOCaug_matrix(1+Np:Np+Nr,1      :Np      ));
  JOC.rr=sparse(JOCaug_matrix(1+Np:Np+Nr,1+Np   :Np+Nr   ));
  JOC.rs=sparse(JOCaug_matrix(1+Np:Np+Nr,1+Np+Nr:Np+Nr+Nr));
  
  
  JOC.sp=sparse(JOCaug_matrix(1+Np+Nr:Np+Nr+Nr,1      :Np));
  JOC.sr=sparse(JOCaug_matrix(1+Np+Nr:Np+Nr+Nr,1+Np   :Np+Nr));
  JOC.ss=sparse(JOCaug_matrix(1+Np+Nr:Np+Nr+Nr,1+Np+Nr:Np+Nr+Nr));

  % size(JOC.pp)
  % size(JOC.pr)
  % size(JOC.ps)

  % size(JOC.rp)
  % size(JOC.rr)
  % size(JOC.rs)

  
  % size(JOC.sp)
  % size(JOC.sr)
  % size(JOC.ss)


  
  
  
  
