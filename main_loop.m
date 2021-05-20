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

verb = 1; % verbosity level: {0,1,2}


grounded_node=1

plot=0


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
  for h_i = 1:4
	       % INCRESING TIME STEP
	       % set here 1:5
    for i=1:4
      N=4*(2^(i-1))
      Nt = N+1;

      str_test=sprintf('%s_h%d_rec%d_N%0.5d_',test_case,h_i,rec,N)
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
      
      %set here [9,10,11]
      for sol=[10];
	% Mesh structure:
	% nodes -> array of nodes coordinates [x y]
	% cells -> array of cells nodes [#nodes node1 node2 node3 ...]
	solver_approach=sol;
	

	ctrl_inner11=ctrl_solver;
	ctrl_inner22=ctrl_solver;
	ctrl_outer=ctrl_solver;

	if (sol==8)
	  ctrl_outer.init('stationary_iterative_methods',1e-5,1000);

	  solvers={'agmg'  ,'direct' ,'krylov'  ,'incomplete'};
	  iters  ={10      ,1        ,10        ,  1          };
	  labels  ={'agmg10','direct','krylov10','incomplete'};

	  % set here from solvers
	  for i=1:length(solvers)
	    ctrl_inner11.init(solvers{i},1e-12,iters{i},1.0,0,labels{i});

	    
	    extras={'block_diag','block_triangular'}
	    for i=1:2
	      controls = struct('indc',grounded_node,...
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
	  ctrl_outer.init('fgmres',1e-5,1000);

	  solvers={'agmg' ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	  iters  ={1      ,10      ,1       ,1        ,10        ,  1          };
	  labels  ={'agmg1','agmg10','direct','krylov10','krylov10','incomplete'};

	  % set here from solvers
	  for i=1:length(solvers)
	    ctrl_inner11.init(solvers{i},1e-12,iters{i},1.0,0,labels{i});
	    extra_info=''; 
	    controls = struct('indc',grounded_node,...
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
	  ctrl_outer.init('fgmres',1e-5,1000);

	  solvers={'direct','agmg' ,'agmg'  ,'incomplete','krylov' ,'krylov'  };
	  iters  ={1       ,100      ,10       ,1           ,100        ,  1          };
	  label  ={'direct','agmg100','agmg10','incomplete','krylov10','krylov10'};

	  % set here from solvers
	  for i=[1,3];%length(solvers)
	    ctrl_inner11.init(solvers{i},1e-13,iters{i},1.0,0,label{i});
	    

	    % cycle approach for S inversion
	    invS_approach={'full','block_triang', 'block_diag'};

	    % set here 
	    for j=1:1%length(invS_approach)
	      controls = struct('indc',grounded_node,...
				'sol',solver_approach,...
				'ctrl_inner11',ctrl_inner11,...
				'ctrl_inner22',ctrl_inner22,...
				'ctrl_outer',ctrl_outer,...
				'compute_eigen',compute_eigen,...
				'verbose',verbose,...
				'extra_info',invS_approach{j});

	      approach_string=def_approach_string(controls);
	      
	      %controls_string=strcat('_',label{i});

	      geod;
	    end
	  end
	elseif (sol==11)
	  % set here bicgstab,gmres,fgmres (for non stationary prec)
	  ctrl_outer.init('bicgstab',1e-5,1000);

	  % set here other approximate inverse of block11
	  ctrl_inner11.init('diag',1e-12,1,1.0,0,'diag');

	  
	  solvers={'agmg' ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
	  iters  ={1      ,10      ,1       ,1        ,10        ,  1          };
	  label  ={'agmg1','agmg10','direct','krylov1','krylov10','incomplete'};

	  % set here from solvers
	  for i=[2];%1:length(solvers)
	    ctrl_inner22.init(solvers{i},1e-12,iters{i},1.0,0,label{i});
	    extra_info='';
	    controls = struct('indc',grounded_node,...
			      'sol',solver_approach,...
			      'ctrl_inner11',ctrl_inner11,...
			      'ctrl_inner22',ctrl_inner22,...
			      'ctrl_outer',ctrl_outer,...
			      'compute_eigen',compute_eigen,...
			      'verbose',verbose,...
			      'extra_info',extra_info);
	    
	    approach_string=def_approach_string(controls);
	    geod;
	  end
	end 
      end
    end
  end
end
