function dl = get_dl(JF,F,d,W_mat) 
  % dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)

  ncellrho=JF.ncellrho;
  N = JF.ntimestep;

  %
  dl=-W_mat*([JF.rp JF.rr JF.rs]*d + F.r);
  for k = 1:N
		dl(k) = dl(k)/norm(W_mat(k,:))^2;
  end
    
  
  if (0)
    for i = 1:N    
      fprintf('dlambda(%d)=%1.2e\n',i,dl(i)) 
    end
  end
end
