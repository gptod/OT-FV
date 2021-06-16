if (0)

    A = sparse(JF.pp); B1T = sparse(JF.pr); B2 = sparse(JF.rp);
    C = sparse(JF.rr - JF.rs*(JF.ss\JF.sr));
    f1 = f;
    f2 = g-JF.rs*(JF.ss\h);

	      % swap C sign for having standard saddle point notation 
    C = -C;
    if (swap_sign)
      A=-A;
      B1T=-B1T;
      B2=-B2;
      C=-C;
      f1=-f1;
      f2=-f2;
    end

				% assembly full system
    jacobian = [A B1T; B2 -C];

				% assembly rhs
    rhs=[f1;f2];

				% grounding
    if (indc>0)
		       % find nodes at time time step diag(A^i) is max
      disp('GROUNDING')
      indeces=set_grounding_node(A,ncellphi);
      inode=indeces(1);
      [A,B1T,rhs]=grounding(A,B1T,rhs,inode,0);
    else
      irow=indeces_global(1)
      A(irow,1:ncellphi) = vectors_x(:,1)';
      B1T(irow,:)=sparse(1,Nr);
      rhs(irow,:)=0.0;
    end

    [A,B1T,rhs]= manipulate_AB1Trhs(A,B1T,rhs,indeces_global,vectors_x,vectors_y,alphas);

    S = A+B1T*(invdiagC*B2);  
    fp = rhs(1:Np)+B1T*(invC.apply(rhs(Np+1:Np+Nr)));
    
    fprintf(' res |S x -f_r|/f_r = %1.4e \n',norm(S*d(1:Np)-fp)/norm(fp))

    jacobian=[A B1T;B2 -C];


    
    fprintf(' res modified = %1.4e \n',norm(jacobian*d-rhs)/norm(rhs))
  end
