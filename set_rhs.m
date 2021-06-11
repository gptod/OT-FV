function rhs_new = set_rhs(rhs,indeces_local)
  rhs_new=rhs;
  Nt=size(indeces_local)
  for i = 1:Nt
    rhs_new(indeces_global(i))=0.0;
  end
end
