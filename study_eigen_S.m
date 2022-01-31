function [Mtt,Mss,dense_S]=study_eigen_S(A,invA,B2_time,B2_space,B1T_time,B1T_space,B1T_perturbation,C,diagrho,LT)
				% A^{-1}B1T=A^{-1} (Dt^T+div L^T+P)
  Np=size(A,1);
  Nr=size(C,1);
  disp('invA_B1T_time')
  invA_B1T_time =prec_times_matrix(invA,B1T_time);
  disp('invA_B1T_space')
  invA_B1T_space=prec_times_matrix(invA,B1T_space);
  disp('invA_B1T_pert')
  invA_B1T_pert =prec_times_matrix(invA,B1T_perturbation);
 
  
  Mtt=(B2_time )*invA_B1T_time;
  Mts=(B2_time )*invA_B1T_space;
  Mst=(B2_space)*invA_B1T_time;
  Mss=(B2_space)*invA_B1T_space;
  Mtp=(B2_time )*invA_B1T_pert;
  Msp=(B2_space)*invA_B1T_pert;


  dense_S=C+Mtt+Mts+Mst+Mss+Mtp+Msp;

  
   %only_space=B2_invA_B1T-only_time-perturbation;
  ctrl_loc=ctrl_solver;
  ctrl_loc.init('direct',1e-14,0,0.0,0,'invS');
  real_inverseS=sparse_inverse;
  real_inverseS.init(dense_S,ctrl_loc);

  if (0)
  M_C           = prec_times_matrix(@(x) real_inverseS.apply(x), C);
  M_time_time   = prec_times_matrix(@(x) real_inverseS.apply(x), Mtt);
  M_space_space = prec_times_matrix(@(x) real_inverseS.apply(x), Mss);
  M_time_space  = prec_times_matrix(@(x) real_inverseS.apply(x), Mts);
  M_space_time  = prec_times_matrix(@(x) real_inverseS.apply(x), Mst);
  M_P           = prec_times_matrix(@(x) real_inverseS.apply(x), Msp+Mtp);

    
  mc=norm(M_C);
  mtt=norm(M_time_time);
  mts=norm(M_time_space);
  mst=norm(M_space_time);
  mss=norm(M_space_space);
  mp=norm(M_P);
  total=M_C+M_time_time+M_time_space+M_space_time+M_space_space+M_P ;
  disp('Norm that should be zero')
  norm(total-speye(Nr,Nr))

  fprintf('mc=%1.1e mtt=%1.1e mts=%1.1e mst=%1.1e mss=%1.1e mp=%1.1e \n',mc,mtt,mts,mst,mss,mp);

  %candidate=M_C+M_time_time+M_time_space+M_time_space+M_space_time+M_P ;
  end 

  disp('candidate')
  prec=C+Mtt+Mss;
  candidate=prec_times_matrix(@(x) real_inverseS.apply(x),prec);
  
  
  [real,img]=study_eigenvalues(prec\dense_S, 'S^{-1}(Candidate)',1);

  figure
  plot([real,img],'o')
  
  
  %figure
  %spy(candidate)
  %saveas(gcf,strcat('candidate.png'));
end 

  
    
