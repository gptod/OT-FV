function []=test_vectors(d,vectors_x,vector_y,alphas,N,Np,Nr)
  ncellphi=Np/(N+1);
  fprintf('TESTING p1^T x_i + p2^t x_{i,i+1} + q^y=alphas(i)\n')
  for i=1:N
    %p1^T x_i + p2^t x_{i,i+1} + q^y=alphas(i)
    v1=vectors_x(:,1+(i-1)*2);
    v2=vectors_x(:,2+(i-1)*2);
    vy=vector_y(:,i);
    res=abs(...
	     v1'*d(1+(i-1)*ncellphi:    i*ncellphi)+...
	     v2'*d(1+i    *ncellphi:(i+1)*ncellphi)+...
	     vy'*d(Np+1:Np+Nr)-...
	     alphas(i) );
    fprintf(' %d alpha=%1.2e  res=%1.2e relres=%1.2e \n', i, alphas(i),res,res/alphas(i))
  end
end  
