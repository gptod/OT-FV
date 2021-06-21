% Time discretization. Uniform grid with Nt=N+1 intervals. N is the number
% of intermediate densities. Choose N odd, N>=1, in order to have the
% approximate geodesic midpoint.
%N = 16;
%test_case='gauss'
clear;
test_case='sin'
				%test_case='sin_old'
				%test_case='cross'

	% Type of reconstruction for the density on the diamond cells.
	% 1 -> weighted arithmetic mean
	% 2 -> weighted harmonic mean
rec = 1;

eps_0 = 1e-6; % tolerance 

verb = 1; % verbosity level: {0,1,2}



plot_figures=0

restart=0;

save_data=1;
read_from_file=0;
h5_file2read='';


compute_eigen=0;
verbose=1;

     
    % Space discretization. Three types of mesh families available:
    % 1 -> regular triangulation of the domain, with only acute angles
%      (https://www.i2m.univ-amu.fr/fvca5/benchmark/Meshes/index.html)
% 2 -> nested-constructed meshes, based on the previous family
%      (find details at https://hal.archives-ouvertes.fr/hal-03032446)
% 3 -> cartesian grids

% set here 1 to 2
for mesh_type = 2
  % INCRESINg DISCRETIZATION SPACE
  % For each mesh, five levels of refinement h_i, 1->5, are available.

  % set here for 1 to 4
  for h_i = 1
	       % INCRESING TIME STEP
	       % set here 1:5
    for i=1
      N=4*(2^(i-1))
      Nt = N+1;

      str_test=sprintf('%s_h%d_rec%d_N%0.5d_',test_case,h_i,rec,N)
      %
      % 9 : P=(A B1T)
      %       (0 -C )

      % 10 : P=(~S B1T)
      %        (0  -C )
      % ~S=A+B1T C^{-1} B2  
      % ~S=upper_triang( A+B1T C^{-1} B2 )
      % ~S=block_diag(   A+B1T C^{-1} B2 )

      % 11 : P=(~A   B1T)
      %        (B2  -~SCA )
      % ~A=diag(A)
      % ~SCA=C+B2 diag(A)^{-1} B1T


      % 12 :  P=(~A   B1T)
      %        (B2  -~SCA )
      %  ~A^{-1} = inverse of A
      %  ~SCA=C+B2 diag(A)^{-1} B1T
      %  (~SCA)^{-1}= approx. inverse
      
      %set here [9,10,11]
      for sol=[12];
	% Mesh structure:
	% nodes -> array of nodes coordinates [x y]
	% cells -> array of cells nodes [#nodes node1 node2 node3 ...]
	solver_approach=sol;
	

	ctrl_inner11=ctrl_solver;
	ctrl_inner22=ctrl_solver;
	ctrl_outer=ctrl_solver;

	if (sol==8)
	  ctrl_outer.init('stationary_iterative_methods',1e-9,1000);

	  solvers={'agmg'  ,'direct' ,'krylov'  ,'incomplete'};
	  iters  ={10      ,1        ,10        ,  1          };
	  labels  ={'agmg10','direct','krylov10','incomplete'};

	  % set here from solvers
	  for i=1:length(solvers)
	    ctrl_inner11.init(solvers{i},1e-12,iters{i},1.0,0,labels{i});

	    
	    extras={'block_diag','block_triangular'}
	    for i=1:2
	      controls = struct('save_data',save_data,...
				'indc',grounded_node,...
				'sol',solver_approach,...
				'ctrl_inner11',ctrl_inner11,...
				'ctrl_inner22',ctrl_inner22,...
				'ctrl_outer',ctrl_outer,...
				'compute_eigen',compute_eigen,...
				'verbose',verbose,...
				'extra_info',extras{i});
	      approach_string=def_approach_string(controls);
	      geod
	    end
	  end
	    
	

	elseif (sol==9)
	    % set here bicgstab,gmres,fgmres (for non stationary prec)
	  ctrl_outer.init('bicgstab',1e-5,100);

	  solvers={'agmg' ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	  iters  ={1      ,10      ,1       ,1        ,10        ,  1          };
	  labels  ={'agmg1','agmg10','direct','krylov10','krylov10','incomplete'};

	  % set here from solvers
	  for i=1:length(solvers)
	    ctrl_inner11.init(solvers{i},1e-12,iters{i},1.0,0,labels{i});
	    extra_info=''; 
	    controls = struct('save_data',save_data,...
			      'indc',grounded_node,...
			      'sol',solver_approach,...
			      'ctrl_inner11',ctrl_inner11,...
			      'ctrl_inner22',ctrl_inner22,...
			      'ctrl_outer',ctrl_outer,...
			      'compute_eigen',compute_eigen,...
			      'verbose',verbose,...
			      'extra_info',extra_info);
	    
	    approach_string=def_approach_string(controls);
	    geod
	  end
	elseif (sol==10)
	  % set here bicgstab,gmres,fgmres (for non stationary prec)
	  ctrl_outer.init('fgmres',1e-5,200,0.0,0);

	  % preconditioner approach
	  outer_prec='full'

          % set here other approximate inverse of block22
	  relax4inv22=0;
	  ctrl_inner22.init('diag',... %approach
			    1e-12,... %tolerance
			    1,...% itermax
			    0.0,... %omega
			    0); %verbose
	  

	  % cycle approach for S inversion
	  % full:        :~S=S
	  % block_triang :~S=upper block triangular of S
	  % block_diag   :~S=block diagonal of S
	  invS_approach={'full','block_triang', 'block_diag'};
	  grounded_node=1;
	  relax4inv11=0;

	  
	  % set here 
	  for j=[1]%length(invS_approach)

	    % set solver for block 11 (schurAC)
	    solvers={'direct','agmg'   ,'agmg'  ,'agmg' ,'incomplete','krylov'  ,'krylov'  };
	    iters  ={1       ,100      ,10      ,1      ,1           ,100       ,1         };
	    label  ={'direct','agmg100','agmg10','agmg1','incomplete','krylov10','krylov10'};

   	    % set here from solvers
	    for i=[4];%length(solvers)
	      ctrl_inner11.init(solvers{i},1e-12,iters{i},1.0,0,label{i});
	      
	      
	      % store all controls
	      controls = struct('save_data',save_data,...
				'indc',grounded_node,...
				'sol',solver_approach,...
				'outer_prec',outer_prec,...
				'ctrl_inner11',ctrl_inner11,...
				'ctrl_inner22',ctrl_inner22,...
				'ctrl_outer',ctrl_outer,...
				'compute_eigen',compute_eigen,...
				'verbose',verbose,...
				'extra_info',invS_approach{j},...
				'relax4inv11',relax4inv11,...
				'relax4inv22',relax4inv22);
	      

	      approach_string=strcat('testprepostprocess_g_schurACwithdiagC_',...
				   ctrl_outer.approach,'_',...
				   outer_prec,'_',...
				   'invSAC',invS_approach{j},ctrl_inner11.label,'_',...
				   'invC',ctrl_inner22.label);
	      

	      if (restart)
		approach_string=strcat('restarted_',approach_string);
	      end

	      disp(approach_string)
              %approach_string=def_approach_string(controls);

	      geod;
	    end
	  end
	elseif (sol==11)
	  % krylov based solvers for M [x;y] = [f;g]
	  % pcg works only with full
	  outer_solvers={'bicgstab'  ,'gmres','fgmres' ,'pcg'};

	  % set here fgmres (for non stationary prec), bicgstab,gmres, pcg
	  for isolver=[3]%1:length(outer_solvers)
	  ctrl_outer.init(outer_solvers{isolver},1e-5,1000,0.0,0); % verbose=[1,2] works only for fgmres

	  
	  % external prec appraoch
	  outer_precs={'full' ,'lower_triang'  ,'upper_triang','identity'};
	  nouter_precs=length(outer_precs);

	  % pcg works only with full
	  if (strcmp(outer_solvers,'pcg'))
	    nouter_precs=1
	  end

	  % fgmres works only with right prec
	  %if (strcmp(outer_solvers,'fgmres'))
	  %  nouter_precs=1
	  %end
	  left_right='right';

	  
	  for iprec=[1,3]%1:nouter_precs
	    outer_prec=outer_precs{iprec};
	  
	  % set here other approximate inverse of block11
	  ctrl_inner11.init('diag',... %approach
			    1e-12,... %tolerance
			    10,...% itermax
			    0.0,... %omega
			    0,... %verbose
			    'diag'); %label
	  %extra_info='full';
	  extra_info='block';
	  relax4_inv11=0e-12;
	  
	  % set grounded_node>0 to gorund the potential in grounded node
	  grounded_node=0;

	  diagonal_scaling=1;

	  % set here list of solvers for block 22 
	  solvers={'agmg' ,'agmg'  ,'agmg' ,'direct','krylov' ,'krylov'  ,'incomplete','diag'};
	  iters  ={1      ,10      ,100,1       ,1        ,10        ,  1          ,0};
	  label  ={'agmg1','agmg10','agmg100','direct','krylov1','krylov10','incomplete','diag'};
	  relax4_inv22=0;
	  
	  for i=[1];%1:length(solvers)
	    ctrl_inner22.init(solvers{i},1e-13,iters{i},1.0,0,label{i});
	    controls = struct('save_data',save_data,...
			      'indc',grounded_node,...
			      'diagonal_scaling',diagonal_scaling,...
			      'sol',solver_approach,...
			      'outer_prec',outer_prec,...
			      'left_right',left_right,...
			      'ctrl_inner11',ctrl_inner11,...
			      'ctrl_inner22',ctrl_inner22,...
			      'ctrl_outer',ctrl_outer,...
			      'compute_eigen',compute_eigen,...
			      'verbose',verbose,...
			      'extra_info',extra_info,...
			      'relax4inv11',relax4_inv11,...
			      'relax4inv22',relax4_inv22);

	    approach_string=strcat('original_schurCAwithdiagA_',...
				   ctrl_outer.approach,'_',...
				   left_right,'_',outer_prec,'_prec_',...
				   'invA',ctrl_inner11.label,'_',...
				   'invSCA',ctrl_inner22.label);
				   %...'_fullprecSCA_invA',label{i},'_invSCA',ctrl_inner22.label);

	    disp(approach_string)
	    %approach_string=def_approach_string(controls);
	    geod;
	  end
	  end
	  end
	elseif (sol==12)
	  % set here bicgstab,gmres,fgmres (for non stationary prec)
	  ctrl_outer.init('fgmres',1e-5,1000,0.0,0);

	  % we can use three approach 
	  outer_precs={'full' ,'lower_triang'  ,'upper_triang','identity'};

	  diagonal_scaling=0;
	  
          % we need to ground the solution since A_11 is singular
	  % grounded<0 C^T x1 =0
	  % grounded=0 no grounding
	  % grounded>0 solution is grounded on one node
	  grounded_node=1;

          % left or right preconditoner
	  % only right for fgmres
	  left_right='right';

	  % we may add an identity matrix to the blocks
	  relax4_inv11=0;
	  relax4_inv22=0;

	  % 
	  compute_eigen=1;
	  
	  for j=[1]
	    outer_prec=outer_precs{j};

	    approach_inverse_A='full';
	    approach_inverse_A='block';

	    
	    solvers11={'diag','agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	    iters11  ={1      ,1      ,10       ,1       ,1        ,100        ,  1         };
	    label11  ={'diagA','agmgA1','agmgA10' ,'directA','krylovA1','krylovA100','incompleteA'};
	    for m=[3]
	      

		       % set here other approximate inverse of block11
	      ctrl_inner11.init(solvers11{m},1e-12,iters11{m},1.0,0,label11{m});

	      approaches_schurCA={'diagA','iterative','full'};
	      for k=[1]%:length(approaches_schurCA)

		approach_schurCA=approaches_schurCA{k};

		strcmp(approach_schurCA,'iterative')
		if (strcmp(approach_schurCA,'diagA'))
		  solvers={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
		  iters  ={4      ,1       ,1       ,1        ,100        ,  1         };
		  label  ={'agmg4','agmg1' ,'direct','krylov1','krylov100','incomplete'};
		elseif (strcmp(approach_schurCA,'full'))
		  solvers={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
		  iters  ={10      ,1       ,1       ,1        ,100        ,  1         };
		  label  ={'agmg10','agmg1' ,'direct','krylov1','krylov100','incomplete'};
		elseif (strcmp(approach_schurCA,'iterative'))
		  solvers={'gmres'  ,'bicgstab'  ,'fgmres','pcg' }
		  iters  ={200      ,100       ,100       , 100 };
		  label  ={'gmres10','bicgstab100' ,'fgmres100','pcg100'};
		end

		
				% set here from solvers
		for i=[1];%1:length(solvers)
		  ctrl_inner22.init(solvers{i},1e-13,iters{i},1.0,0,strcat(label{i},'_S'));
		  extra_info='block';
		  controls = struct('save_data',save_data,...
				    'indc',grounded_node,...
				    'sol',solver_approach,...
				    'outer_prec',outer_prec,...
				    'left_right',left_right,...
				    'diagonal_scaling',diagonal_scaling,...
				    'approach_schurCA',approach_schurCA,...
				    'approach_inverse_A',approach_inverse_A,...
				    'ctrl_inner11',ctrl_inner11,...
				    'ctrl_inner22',ctrl_inner22,...
				    'ctrl_outer',ctrl_outer,...
				    'compute_eigen',compute_eigen,...
				    'verbose',verbose,...
				    'extra_info',extra_info,...
				    'relax4inv11',relax4_inv11,...
				    'relax4inv22',relax4_inv22);

		  approach_string=strcat(extra_info,'_manipulated_',ctrl_outer.approach,...
					 'invA',ctrl_inner11.approach,'_precSCA',label{i},...
					 '_invSCA_',approach_schurCA,'_',ctrl_inner22.label);

		  geod;
		end
	      end
		end
	  end
	elseif (sol==13)
	  % set here bicgstab,gmres,fgmres (for non stationary prec)
	  ctrl_outer.init('fgmres',1e-5,300);

          % left or right prec.
	  left_right='right'
	  
          % handle singularity
	  % change A,B1T,f to remove part of its kernel
	  manipulate=0
	  % we need to ground the solution since A_11 is singular
	  % grounded<0 C^T x1 =0
	  % grounded=0 no grounding
	  % grounded>0 solution is grounded on one node
	  grounded_node=1;

          % alpha relaxation
	  alpha=0.1;

				%
	  ctrl_innerA=ctrl_solver;
	  ctrl_innerC=ctrl_solver;
	  ctrl_innerS=ctrl_solver;

	  approach_inverse_A='block';
	  approach_inverse_A='full';

	  approach_inverse_S='SCA';
	  %approach_inverse_S='SAC';

	  approach_prec='SH';
	  %approach_prec='HS';

	  diagonal_scaling=1;
	  
	  
	  solversA={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	  itersA  ={1      ,10       ,1       ,1        ,100        ,  1         };
	  labelA  ={'agmg1','agmg1' ,'direct','krylov1','krylov100','incomplete'};

	  solversS={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	  itersS  ={1      ,10       ,1       ,1        ,100        ,  1         };
	  labelS  ={'agmg1','agmg1' ,'direct','krylov1','krylov100','incomplete'};

	  ctrl_innerC.init('diag',1e-12,1,1.0,0,'C');
	  
          % IMPORTANT: If we want to use 2 agmg solvers we have to set preprocess=0
	  % in sparse inverse 
	  for i=[1];%1:length(solversA)
		    % set here other approximate inverse of A
	    
	    ctrl_innerA.init(solversA{i},1e-12,itersA{i},1.0,0,labelA{i});

	    for j=[1];%1:length(solversA)
		      % set here other approximate inverse of A
	      
	      ctrl_innerS.init(solversS{j},1e-12,itersS{j},1.0,0,labelS{j});


	      
	      controls = struct('save_data',save_data,...
				'indc',grounded_node,...
				'sol',solver_approach,...
				'alpha',alpha,...
				'approach_prec',approach_prec,...
				'manipulate',manipulate,...
				'ctrl_outer',ctrl_outer,...
				'left_right',left_right,...
				'diagonal_scaling',diagonal_scaling,...
				'approach_inverse_S',approach_inverse_S,...
				'approach_inverse_A',approach_inverse_A,...
				'ctrl_innerA',ctrl_innerA,...
				'ctrl_innerC',ctrl_innerC,...
				'ctrl_innerS',ctrl_innerS,...
				'compute_eigen',compute_eigen,...
				'verbose',verbose);

	      approach_string=strcat(approach_prec,'_invA',labelA{i},'_invS_',approach_inverse_S,labelS{j});
	      disp(approach_string)
	      geod;	    
	    end
	  end
	end
      end
    end
  end
end
