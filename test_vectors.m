for i=1:N
  v1=vectors(:,1+(i-1)*2);
  v2=vectors(:,2+(i-1)*2);
  res=abs(...
	   v1'*d(1+(i-1)*ncellphi:    i*ncellphi)+...
	   v2'*d(1+i    *ncellphi:(i+1)*ncellphi)-alphas(i) );
  fprintf(' %d alpha=%f  res=%f\n', i, alphas(i),res)
end
  
