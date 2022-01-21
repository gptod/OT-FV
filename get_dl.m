function dl = get_dl(JF,F,d) 
  % dl=-(deltat*|m|)^{-2} W_mat*(B*x + R*y + M*z+g)

  ncellrho=JF.ncellrho;
  N = JF.ntimestep;
  area2h=JF.area2h;

  %
  res_r=-([JF.rp JF.rr JF.rs]*d + F.r);
  factor=1/(norm(area2h)^2);
  
  dl=zeros(N,1);
  for i = 1:N    
    dl(i) = factor*(area2h)'*res_r(1+(i-1)*ncellrho:i*ncellrho);
  end
    
  
  if (0)
    for i = 1:N    
      fprintf('dlambda(%d)=%1.2e\n',i,dl(i)) 
    end
  end
end
