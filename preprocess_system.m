function [A,B1T,rhs,B1T_perturbation]=preprocess_system(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas,controls)
  indc=controls.indc;
  Nr=size(B1T,2);
  N=size(indeces_global,1)-1;
  ncellphi=size(A,1)/(N+1);

  print=0;
  B1T_perturbation=sparse(size(B1T,1),size(B1T,2));
  if (controls.manipulate > 0 )
    if (print)
      disp('MANIPULATE')
    end
    [A,B1T,rhs,B1T_perturbation]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas);
  end

  
  if (floor(indc) == indc) 
    if (indc>0)
      if (print)
	disp('GROUNDING')
      end
      irow=indeces_global(1);
      [A,B1T,rhs]=grounding(A,B1T,rhs,irow,0);
      B1T_perturbation(irow,:)=sparse(1,Nr);
    elseif(indc<0)
      if (print)
	disp('FIX SUM A11')
      end
      irow=indeces_global(1);
      A(irow,1:ncellphi) = ones(1,ncellphi)/ncellphi;
      B1T(irow,:)=sparse(1,Nr);
      rhs(irow,:)=0.0;
      B1T_perturbation(irow,:)=sparse(1,Nr);
    end
  else
    if (print)
      disp('ADDING eps to A11')
    end
    A(1:ncellphi,1:ncellphi)=  A(1:ncellphi,1:ncellphi)+indc*speye(ncellphi,ncellphi);
  end
  
 

  
