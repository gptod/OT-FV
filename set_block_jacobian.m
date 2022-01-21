function [JOC] = set_block_jacobian(JOC,variable_row, variable_col, irow, icol, vector, line,add)
  % JOC :structure with the jacobian
  % variable row (phi,rho,s)
  % variable col (phi,rho,s)
  % irow: number of row block
  % icol: number of col block
  % line: local line number were store the vector
  % vector: vector to be included
  % add: 0 no, 1 : add to existing jacobian
  
  ncell=JOC.ncellphi;
  ncell2h=JOC.ncellrho;



  if (strcmp(variable_row,'phi'))
    if (strcmp(variable_col,'phi'))
      row=(irow-1)*ncell+line;
      l_col=(icol-1)*ncell+1;
      r_col=(icol)*ncell;
      JOC.pp=set_line(JOC.pp,row,lcol,rcol,vector,add);
    end
    if (strcmp(variable_col,'rho'))
      row=(irow-1)*ncell+line;
      l_col=(icol-1)*ncell2h+1;
      r_col=(icol)*ncell2h;
      JOC.pr=set_line(JOC.pp,row,lcol,rcol,vector,add);
    end
    if (strcmp(variable_col,'s'))
      row=(irow-1)*ncell+line;
      l_col=(icol-1)*ncell2h+1;
      r_col=(icol)*ncell2h;
      JOC.ps=set_line(JOC.ps,row,lcol,rcol,vector,add);
    end
  end

    
  if (strcmp(variable_row,'rho'))
    if (strcmp(variable_col,'phi'))
      row=(irow-1)*ncell2h+line;
      l_col=(icol-1)*ncell+1;
      r_col=(icol)*ncell;
      JOC.rp=set_line(JOC.rp,row,lcol,rcol,vector,add);
    end
    if (strcmp(variable_col,'rho'))
      row=(irow-1)*ncell2h+line;
      l_col=(icol-1)*ncell2h+1;
      r_col=(icol)*ncell2h;
      JOC.rr=set_line(JOC.rr,row,lcol,rcol,vector,add);
    end
    if (strcmp(variable_col,'s'))
      row=(irow-1)*ncell2h+line;
      l_col=(icol-1)*ncell2h+1;
      r_col=(icol)*ncell2h;
      JOC.rs=set_line(JOC.rs,row,lcol,rcol,vector,add);
    end
  end
  
  if (strcmp(variable_row,'s'))
    if (strcmp(variable_col,'phi'))
      row=(irow-1)*ncell2h+line;
      l_col=(icol-1)*ncell+1;
      r_col=(icol)*ncell;
      JOC.sp=set_line(JOC.sp,row,lcol,rcol,vector,add);
    end
    if (strcmp(variable_col,'rho'))
      row=(irow-1)*ncell2h+line;
      l_col=(icol-1)*ncell2h+1;
      r_col=(icol)*ncell2h;
      JOC.sr=set_line(JOC.sr,row,lcol,rcol,vector,add);
    end
    if (strcmp(variable_col,'s'))
      row=(irow-1)*ncell2h+line;
      l_col=(icol-1)*ncell2h+1;
      r_col=(icol)*ncell2h;
      JOC.ss=set_line(JOC.ss,row,lcol,rcol,vector,add);
    end
  end
end


