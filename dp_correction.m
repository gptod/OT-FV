function y = dp_correction(dp,JF,dl) 
	N = JF.ntimestep;
	deltat=1/(N+1);
	ncellphi=JF.ncellphi;
	
	sum_dl=0;
	y=dp;
  for i = 1:N
    sum_dl=sum_dl-deltat*dl(i);
		slot=i+1;
    y(1+(slot-1)*ncellphi:slot*ncellphi)=y(1+(slot-1)*ncellphi:slot*ncellphi)+sum_dl;
  end
end
