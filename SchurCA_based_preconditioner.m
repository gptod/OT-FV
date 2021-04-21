function [out] = SchurCA_based_preconditioner(in, invA,invSCA,B1T,B2)
  % P out = in
  % out = (x;y)
  % in  = (f;g)
  % P   = (A  B1T)
  %       (B2 -C )
  % SCA=-(C+B2 A^{-1} B1)  
  n=size(B1T,1);
  m=size(B1T,2);
  out=zeros(size(in));

  % t=A^{-1} f
  t=invA(in(1:n));
  % v= B2 * A^{-1} f 
  v=B2*t;
  % w= -v+g
  w=-v+in(1+n:n+m);
  % y = SCA^{-1} w
  %   = SCA^{-1} ( -B2 A^{-1} f)
  out(n+1:n+m)=invSCA(w);

  % w= B1T*y
  w      = B1T * out(1+n:n+m);
  % x = t - A^{-1} w
  out(1:n) = invA ( w );
  out(1:n) = t-out(1:n);
  
