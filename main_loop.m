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
recs={1,2};
for kk=1
  
rec = recs{kk};

eps_0 = 1e-6; % tolerance 

verb = 1; % verbosity level: {0,1,2}



plot_figures=0

restart=0;

save_data=1;
read_from_file=0;
%h5_file2read='runs/sol10/PhiRhoSMuThetasin_h1_rec1_N00032__schurACwithdiagC_fgmres_full_invSACfullagmg1e-1_invC1_diag.h5';

compute_eigen=0;
verbose=0;

mkdir 'runs'

     
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
  for h_i =4
	       % INCRESING TIME STEP
	       % set here 1:5
    for i=5:6
      N=4*(2^(i-1))
      Nt = N+1;

      str_test=sprintf('%s_mesh%d_h%d_rec%d_N%0.5d_',test_case,mesh_type,h_i,rec,N)

      h5_file2read=strcat('runs/sol10/PhiRhoSMuTheta',str_test(1:strlength(str_test)-1),'.h5')
      if (read_from_file>0)
	str=sprintf('restart%d_',read_from_file);
	str_test=strcat(str,str_test);
      end
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
      for sol=[10];

	folder_run=sprintf('run/sol%s',sol);
	mkdir folder_run
	% Mesh structure:
	% nodutees -> array of nodes coordinates [x y]
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
	  ctrl_outer.init('fgmres',1e-5,200,0.0,2);

	  % preconditioner approach
	  outer_prec='full'

          % set here other approximate inverse of block22
	  relax4inv22=0;
	 
	  
	  for kk=1
	  inverseC_approach=1;
	  
	  % cycle approach for S inversion
	  % full:        :~S=S
	  % block_triang :~S=upper block triangular of S
	  % block_diag   :~S=block diagonal of S
	  invS_approach={'full','block_triang', 'block_diag'};

	  grounded_node=1;
	  manipulation_approach=1;
	  manipulate=0;

	  
	  relax4inv11=0;

	  
	  % set here 
	  for j=[1]%length(invS_approach)

	    % set solver for block 11 (schurAC)
	    solvers={'direct','agmg'   ,'agmg'  ,'agmg' ,'incomplete','krylov'  ,'krylov'  };
	    iters  ={1       ,20      ,10      ,1      ,1           ,100       ,1         };
	    label  ={'direct','agmg1e-1','agmg10','agmg1','incomplete','krylov10','krylov10'};

   	    % set here from solvers
	    for i=[2];%length(solvers)
	      ctrl_inner11.init(solvers{i},1e-1,iters{i},1.0,1,label{i});

	      ctrl_inner22.init('diag',... %approach
				1e-6,... %tolerance
				100,...% itermax
				0.0,... %omega
				0,'diag'); %verbose
	      
	      % store all controls
	      controls = struct('save_data',save_data,...
				'indc',grounded_node,...
				'manipulate',manipulate,...
				'manipulation_approach',manipulation_approach,...
				'sol',solver_approach,...
				'outer_prec',outer_prec,...
				'inverseC_approach',inverseC_approach,...
				'ctrl_inner11',ctrl_inner11,...
				'ctrl_inner22',ctrl_inner22,...
				'ctrl_outer',ctrl_outer,...
				'compute_eigen',compute_eigen,...
				'verbose',verbose,...
				'extra_info',invS_approach{j},...
				'relax4inv11',relax4inv11,...
				'relax4inv22',relax4inv22);
	      

	      approach_string=strcat('schurACwithdiagC_',...
				   ctrl_outer.approach,'_',...
				   outer_prec,'_',...
				   'invSAC',invS_approach{j},ctrl_inner11.label,'_',...
				   'invC',num2str(inverseC_approach),'_',ctrl_inner22.label);
	      

	      if (restart)
		approach_string=strcat('restarted_',approach_string);
	      end

	      disp(approach_string)
              %approach_string=def_approach_string(controls);

	      geod;
	    end
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

	  
	  for iprec=[1]%1:nouter_precs
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
	  diagonal_scaling=0;
	  manipulate=0;
	  manipulation_approach=1;

	  

	  

	  % set here list of solvers for block 22 
	  solvers={'agmg' ,'agmg'  ,'agmg' ,'direct','krylov' ,'krylov'  ,'incomplete','diag','krylov_no_prec'};
	  iters  ={1      ,10      ,100,1       ,1        ,10        ,  1          ,0,10};
	  label  ={'agmg1','agmg10','agmg100','direct','krylov1','krylov10','incomplete','diag','purekrylov'};
	  relax4_inv22=0;
	  
	  for i=[1];%1:length(solvers)
	    ctrl_inner22.init(solvers{i},1e-5,iters{i},1.0,0,label{i});
	    controls = struct('save_data',save_data,...
			      'indc',grounded_node,...
			      'diagonal_scaling',diagonal_scaling,...
			      'manipulate',manipulate,...
			      'manipulation_approach',manipulation_approach,...
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

	    approach_string=strcat('grounded',num2str(grounded_node),...
				   '_diagscal',num2str(diagonal_scaling),...
				   '_manA',num2str(manipulate),...
				   '_schurCAwithdiagA_',...
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


	  ctrl_inner22inner=ctrl_solver;
	  
	  
          % we need to ground the solution since A_11 is singular
	  % grounded<0 C^T x1 =0
	  % grounded=0 no grounding
	  % grounded>0 solution is grounded on one node
	  for ii=[-1]
	    for jj=[0]
	      for kk=[1]
	  grounded_node=ii;
	  diagonal_scaling=jj;
	  manipulate=1;
	  manipulation_approach=1;

          % left or right preconditoner
	  % only right for fgmres
	  left_right='right';

	  % we may add an identity matrix to the blocks
	  relax4_inv11=0;
	  relax4_inv22=0;

	  % 
	  compute_eigen=1;
	  
	  for j=[3]
	    outer_prec=outer_precs{j};

	    %approach_inverse_A='full';
	    approach_inverse_A ='block';

	    
	    solvers11={'diag','agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	    iters11  ={1      ,1      ,100       ,1       ,1        ,100        ,  1         };
	    label11  ={'diagA','agmgA1','agmg1e-1' ,'directA','krylovA1','krylovA100','incompleteA'};
	    for m=[3]
	      

	      % set here other approximate inverse of block11
	      ctrl_inner11.init(solvers11{m},1e-1,iters11{m},1.0,0,label11{m});

	      approaches_schurCA={'diagA','iterative','full','iterative+SwithdiagA','timelaplacian'};
	      for k=[1]%:length(approaches_schurCA)

		approach_schurCA=approaches_schurCA{k};

		if (strcmp(approach_schurCA,'diagA'))
		  solvers22={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
		  iters22  ={10      ,1       ,1       ,1        ,100        ,  1         };
		  label22  ={'agmg1e-1','agmg1' ,'direct','krylov1','krylov100','incomplete'};
		elseif (strcmp(approach_schurCA,'full'))
		  solvers22={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
		  iters22  ={10      ,1       ,1       ,1        ,100        ,  1         };
		  label22  ={'agmg10','agmg1' ,'direct','krylov1','krylov100','incomplete'};
		elseif (strcmp(approach_schurCA,'iterative'))
		  solvers22={'gmres'  ,'bicgstab'  ,'fgmres','pcg' }
		  iters22  ={200      ,100       ,100       , 100 };
		  label22  ={'gmres10','bicgstab100' ,'fgmres100','pcg100'};
		elseif (strcmp(approach_schurCA,'iterative+SwithdiagA'))
		  solvers22={'gmres'  ,'bicgstab'  ,'fgmres','pcg' }
		  iters22  ={200      ,100       ,100       , 100 };
		  label22  ={'gmres10','bicgstab100' ,'fgmres100','pcg100'};
		  
		  
		  ctrl_inner22inner.init('agmg',1e-2,100,1.0,0,'SCAwithdiagA');
		elseif (strcmp(approach_schurCA,'timelaplacian'))
		  solvers22={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
		  iters22  ={10      ,1       ,1       ,1        ,100        ,  1         };
		  label22  ={'agmg10','agmg1' ,'direct','krylov1','krylov100','incomplete'};
		end

		
				% set here from solvers
		for i=[1];%1:length(solvers)
		  ctrl_inner22.init(solvers22{i},1e-1,iters22{i},1.0,0,strcat(label22{i},'_S'));
		  extra_info='block';
		  controls = struct('save_data',save_data,...
				    'indc',grounded_node,...
				    'sol',solver_approach,...
				    'outer_prec',outer_prec,...
				    'left_right',left_right,...
				    'diagonal_scaling',diagonal_scaling,...
				    'manipulate',manipulate,...
				    'manipulation_approach',manipulation_approach,...
				    'approach_schurCA',approach_schurCA,...
				    'approach_inverse_A',approach_inverse_A,...
				    'ctrl_inner11',ctrl_inner11,...
				    'ctrl_inner22',ctrl_inner22,...
				    'ctrl_inner22inner',ctrl_inner22inner,...
				    'ctrl_outer',ctrl_outer,...
				    'compute_eigen',compute_eigen,...
				    'verbose',verbose,...
				    'extra_info',extra_info,...
				    'relax4inv11',relax4_inv11,...
				    'relax4inv22',relax4_inv22);

		  approach_string=strcat('grounded',num2str(grounded_node),...
					 '_diagscal',num2str(diagonal_scaling),...
					 '_manA',num2str(manipulate),'_',...
					 'sol12_',...
					 ctrl_outer.approach,'_',...
					 outer_prec,'_',...
					 'invA_',approach_inverse_A,label11{m},...
					 '_precSCA_',approach_schurCA,'_',label22{i})

		  geod;
		end
	      end
	    end

	  end
	      end
	    end
	    
	  end
	elseif (sol==13)
	  % set here bicgstab,gmres,fgmres (for non stationary prec)
	  ctrl_outer.init('fgmres',1e-5,1000);

          % left or right prec.
	  left_right='right'
	  
          % handle singularity
	  % change A,B1T,f to remove part of its kernel
	  manipulate=0;
	  manipulation_approach=1;
	  % we need to ground the solution since A_11 is singular
	  % grounded<0 C^T x1 =0
	  % grounded=0 no grounding
	  % grounded>0 solution is grounded on one node
	  grounded_node=0;
	  diagonal_scaling=1;


				% alpha relaxation
	  alphas={0.01,0.05,0.1,0.5,0.9};
	  for ll=4
	  alpha=alphas{ll};

				%
	  ctrl_innerA=ctrl_solver;
	  ctrl_innerC=ctrl_solver;
	  ctrl_innerS=ctrl_solver;

	  %approach_inverse_A='block';
	  approach_inverse_A='full';

	  approach_inverse_S='SCA';
	  %approach_inverse_S='SAC';

	  approach_prec='SH';
	  %approach_prec='HS';

	  	  
	  
	  solversA={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	  itersA  ={1      ,100       ,1       ,1        ,100        ,  1         };
	  labelA  ={'agmg1','agmg1e-1' ,'direct','krylov1','krylov100','incomplete'};

	  solversS={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	  itersS  ={1      ,100       ,1       ,1        ,100        ,  1         };
	  labelS  ={'agmg1','agmg1e-1' ,'direct','krylov1','krylov100','incomplete'};

	  ctrl_innerC.init('agmg',1e-1,100,1.0,0,'C');
	  
          % IMPORTANT: If we want to use 2 agmg solvers we have to set preprocess=0
	  % in sparse inverse 
	  for i=[2];%1:length(solversA)
		    % set here other approximate inverse of A
	    
	    ctrl_innerA.init(solversA{i},1e-1,itersA{i},1.0,0,labelA{i});

	    for j=[2];%1:length(solversA)
		      % set here other approximate inverse of A
	      
	      ctrl_innerS.init(solversS{j},1e-1,itersS{j},1.0,0,labelS{j});


	      
	      controls = struct('save_data',save_data,...
				'indc',grounded_node,...
				'sol',solver_approach,...
				'alpha',alpha,...
				'approach_prec',approach_prec,...
				'manipulate',manipulate,...
				'manipulation_approach',manipulation_approach,...
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

	      approach_string=strcat(approach_prec,'alpha',num2str(alpha),'_invA',labelA{i},'_invS_',approach_inverse_S,labelS{j});
	      disp(approach_string)
	      geod;	    
	    end
	  end
	  end
	elseif (sol==14)
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
	  
	  for kk=1
	  inverseC_approach=kk;
	  
	  % cycle approach for S inversion
	  % full:        :~S=S
	  % block_triang :~S=upper block triangular of S
	  % block_diag   :~S=block diagonal of S
	  invS_approach={'full'};

	  grounded_node=1;
	  manipulation_approach=1;
	  manipulate=0;

	  W_approach='C';
				%W_approach='cutC';
	  inverse22_approach='diag';
	  lower_bound=1e-12;
	  upper_bound=1e12;
	  
	  
	  relax4inv11=0;

	  
	  % set here 
	  for j=[1]%length(invS_approach)

	    % set solver for block 11 (schurAC)
	    solvers={'direct','agmg'   ,'agmg'  ,'agmg' ,'incomplete','krylov'  ,'krylov'  };
	    iters  ={1       ,100      ,10      ,1      ,1           ,100       ,1         };
	    label  ={'direct','agmg1e-1','agmg10','agmg1','incomplete','krylov10','krylov10'};

   	    % set here from solvers
	    for i=[2];%length(solvers)
	      ctrl_inner11.init(solvers{i},1e-1,iters{i},1.0,0,label{i});
	      
	      
	      % store all controls
	      controls = struct('save_data',save_data,...
				'indc',grounded_node,...
				'manipulate',manipulate,...
				'manipulation_approach',manipulation_approach,...
				'sol',solver_approach,...
				'outer_prec',outer_prec,...
				'inverseC_approach',inverseC_approach,...
				'W_approach',W_approach,...
				'lower_bound',lower_bound,...
				'upper_bound',upper_bound,...
				'ctrl_inner11',ctrl_inner11,...
				'ctrl_inner22',ctrl_inner22,...
				'ctrl_outer',ctrl_outer,...
				'compute_eigen',compute_eigen,...
				'verbose',verbose,...
				'extra_info',invS_approach{j},...
				'relax4inv11',relax4inv11,...
				'relax4inv22',relax4inv22);
	      

	      approach_string=strcat('schurACwithdiagC_',...
				   ctrl_outer.approach,'_',...
				   outer_prec,'_',...
				   'invSAC',invS_approach{j},ctrl_inner11.label,'_',...
				   'invC',num2str(inverseC_approach),'_',ctrl_inner22.label);
	      

	      
	      disp(approach_string)
              
	      geod;
	    end
	  end
	  end
	end
	
      end
    end
  end
end
end
