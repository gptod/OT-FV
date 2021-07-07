function [out] = SchurAC_based_preconditioner(in, invSAC,invC,B1T,B2,prec_type,dimblock)
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

  n=size(B1T,1);
  m=size(B1T,2);

  
  out=zeros(size(in));

  if ( strcmp(type,'full'))
    % P   = (A  B1T)
    %       (B2 -C )  
    % SAC=A+B1T*C^{-1}B2

    % P^{-1} = ( I         0 ) ( SAC^{-1}   0      ) ( I B1^{T} C^{-1} )
    %          ( C^{-1} B2 I ) ( 0          -C^{-1}) ( 0 I             )
   
    % t=C^{-1} g %store this vector
    t=invC(in(1+n:n+m));

    % x = SAC^{-1} ( f + B1T C^{-1} g)
    v= in(1:n) + B1T*t;
    if( verbose)
      disp(n/dimblock)
      for i=1:n/dimblock
	fprintf('imb=%1.4e\n',sum(v(1+(i-1)*dimblock:i*dimblock)))
      end
    end
    out(1:n)=invSAC(v);

    % y = C^{-1} B2 x - C^{-1} g
    %   = C^{-1} B2 x - t
    out(n+1:n+m) = invC(B2*out(1:n))-t;

  elseif(strcmp(type,'lower_triang'))
    % P   = (S     )
    %       (B2  -C)  
    % SAC=A+B1T*C^{-1}B2

  
    if( verbose && dimblock>0 )
      disp(n/dimblock);
      for i=1:n/dimblock
	fprintf('imb=%1.4e\n',sum(in(1+(i-1)*dimblock:i*dimblock)))
      end
    end
    
    
    % x=A^{-1} f
    out(1:n)=invSAC(in(1:n));

    % y = -C^{-1} ( -B2 x + g)
    out(n+1:n+m)=-invC( B2*out(1:n) + in(1+n:n+m));
  elseif(strcmp(type,'upper_triang'))
     % P   = (S  B1T)
     %       (   -C  )  
     % SAC=A+B1T*C^{-1}B2

     % y  = -C^{-1} (  g)
     out(n+1:n+m)=-invC(in(n+1:n+m));
    
     % x = SAC^{-1} (f- B1T y)
     v=in(1:n) - B1T*out(n+1:n+m);
     
     if( verbose && dimblock>0 )
       disp(n/dimblock);
       for i=1:n/dimblock
	 fprintf('imb=%1.4e\n',sum(v(1+(i-1)*dimblock:i*dimblock)))
       end
     end
     out(1:n) = invSAC(v);
   
  end 

end
