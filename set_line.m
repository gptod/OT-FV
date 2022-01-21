function [matrix]=set_line(matrix,row,lcol,rcol,vector,add)
  if (not(add))
    matrix(row,lcol:rcol)=sparse(1:matrix(rcol-lcol+1));
  end
  matrix(row,lcol:rcol)=matrix(row,lcol:rcol)+vector;
end
