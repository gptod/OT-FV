function dl = get_dl(JF,F,d) 
  % dl=-1/area E*(B*x + R*y + M*z-g)
  ncellrho=JF.ncellrho;
  N = JF.ntimestep;

  %
	dl=zeros(N,1);	
	res=[JF.rp JF.rr JF.rs]*d + F.r;
  for k = 1:N
		dl(k) = -sum(res((k-1)*ncellrho+1:k*ncellrho))/JF.areadomain;
  end
    
  
  if (0)
    for k = 1:N    
      fprintf('dlambda(%d)=%1.2e\n',k,dl(k)) 
    end
  end
end
