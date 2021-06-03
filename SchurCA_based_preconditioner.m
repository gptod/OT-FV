function [out] = SchurCA_based_preconditioner(in, invA,invSCA,B1T,B2,prec_type,dimblock)
  % P out = in
  % out = (x;y)
  % in  = (f;g)


  verbose=0;
  if (exist('prec_type','var') )
    type=prec_type;
  else
    type='full';
  end

  if (exist('prec_type','var') )
    dimblock=dimblock;
  else
    dimblock=0;
  end

  n=size(B2,2);
  m=size(B2,1);

  %fprintf('norm v in P^{-1} v - %1.4e \n', norm(in))


  
  out=zeros(size(in));

  if ( strcmp(type,'full'))
    % P   = (A  B1T)
    %       (B2 -C )  
    % SCA=-(C+B2 A^{-1} B1)
    
    
    
    % t=A^{-1} f
    if (verbose)
      disp('first')
      print_imbalance(in(1:n),dimblock)
    end
    t=invA(in(1:n));
    % v= B2 * A^{-1} f 
    v=B2*t;
    % w= -v+g
    w=-v+in(1+n:n+m);
    % y = SCA^{-1} w
    %   = SCA^{-1} ( -B2 A^{-1} f)
    out(n+1:n+m)=invSCA(w);

    % w= B1T*y
    w      = B1T ( out(1+n:n+m) );
    % x = t - A^{-1} w
    if (verbose)
      disp('second')
      print_imbalance(w,dimblock)
    end
    out(1:n) = invA ( w );
    out(1:n) = t-out(1:n);

  elseif(strcmp(type,'lower_triang'))
     % P   = (A     )
     %       (B2  S )  
     % SCA=-(C+B2 A^{-1} B1)
    
     %x=A^{-1} f
     if (verbose)
       disp('first')
       print_imbalance(in(1:n),dimblock)
     end
     out(1:n)=invA(in(1:n));

     
     % w = g - B2 * x 
     w   = in(1+n:n+m) - B2*out(1:n);
     % y = SCA^{-1} w
     %   = SCA^{-1} ( -B2 A^{-1} f+g)
     out(n+1:n+m)=invSCA(w);
  elseif(strcmp(type,'upper_triang'))
    % P   = (A  B1T)
    %       (   S  )  
    % SCA=-(C+B2 A^{-1} B1)
    
    % y  = SCA^{-1} (  g)
    
    out(n+1:n+m)=invSCA(in(n+1:n+m));
    
    
    % v= f - B1T*y
    %print_imbalance(in(1:n),dimblock)
    v=in(1:n) - B1T ( out(n+1:n+m));
    %if (verbose)
    %print_imbalance(v,dimblock)
    %end
    
    % x = A^{-1} v
    out(1:n) = invA(v);

  elseif(strcmp(type,'identity'))
    out=in;
  end 

end

