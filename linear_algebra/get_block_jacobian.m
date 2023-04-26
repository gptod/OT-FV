function [matrix] = get_block_jacobian(JOC,variable_row, variable_col, irow, icol)

  ncell=JOC.ncellphi;
  ncell2h=JOC.ncellrho;



  if (strcmp(variable_row,'phi'))
    if (strcmp(variable_col,'phi'))
      matrix=JOC.pp(1+(irow-1)*ncell  :(irow)*ncell,...
		    1+(icol-1)*ncell  :(icol)*ncell  );
    end
    if (strcmp(variable_col,'rho'))
      matrix=JOC.pr(1+(irow-1)*ncell    :(irow)*ncell,...
		    1+(icol-1)*ncell2h  :(icol)*ncell2h  );
    end
    if (strcmp(variable_col,'s'))
      matrix=JOC.ps(1+(irow-1)*ncell    :(irow)*ncell,...
		    1+(icol-1)*ncell2h  :(icol)*ncell2h  );
    end
  end

    
  if (strcmp(variable_row,'rho'))
    if (strcmp(variable_col,'phi'))
      matrix=JOC.rp(1+(irow-1)*ncell2h  :(irow)*ncell2h,...
		    1+(icol-1)*ncell    :(icol)*ncell  );
    end
    if (strcmp(variable_col,'rho'))
      matrix=JOC.rr(1+(irow-1)*ncell2h  :(irow)*ncell2h,...
		    1+(icol-1)*ncell2h  :(icol)*ncell2h);
    end
    if (strcmp(variable_col,'s'))
      matrix=JOC.rs(1+(irow-1)*ncell2h  :(irow)*ncell2h,...
		    1+(icol-1)*ncell2h  :(icol)*ncell2h);
    end
  end
  
  if (strcmp(variable_row,'s'))
    if (strcmp(variable_col,'phi'))
      matrix=JOC.sp(1+(irow-1)*ncell2h  :(irow)*ncell2h,...
		    1+(icol-1)*ncell    :(icol)*ncell  );
    end
    if (strcmp(variable_col,'rho'))
      matrix=JOC.sr(1+(irow-1)*ncell2h  :(irow)*ncell2h,...
		    1+(icol-1)*ncell2h  :(icol)*ncell2h  );
    end
    if (strcmp(variable_col,'s'))
      matrix=JOC.ss(1+(irow-1)*ncell2h  :(irow)*ncell2h,...
		    1+(icol-1)*ncell2h  :(icol)*ncell2h  );

    end
  end
end


