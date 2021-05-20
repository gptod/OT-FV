function [str] = def_approach_string(controls)

  sol=controls.sol;
  if ( sol==9 )
    str=strcat('krylov_triangprec_invA',controls.ctrl_inner11.label,'_invC');
  elseif(sol==10)
    if ( strcmp(controls.extra_info,'full'))
      str=strcat('krylov_triangprec_invSAC',controls.ctrl_inner11.label,'_invC');
    elseif ( strcmp(controls.extra_info,'block_triang'))
      %str='krylov_triangprec_invblocktriangSAC_invC';
      str=strcat('krylov_triangprec_invblocktriangSAC',controls.ctrl_inner11.label,'_invC');
    elseif ( strcmp(controls.extra_info,'block_diag'))
      str=strcat('krylov_triangprec_invblockdiagSAC',controls.ctrl_inner11.label,'_invC');
      %str='krylov_triangprec_invblockdiagSAC_invC';
    else
      disp('NOT DEFINED',block_diag) 
    end
  elseif(sol==11)
    str=strcat('krylov_fullprec_invA',controls.ctrl_inner11.label,'_invC+B2invdiagAB1T',controls.ctrl_inner22.label);
    %str='invA_invSCA';
  end


      
