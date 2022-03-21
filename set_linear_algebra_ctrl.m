function [ctrls, approach_descriptions] = ...
				 set_linear_algebra_ctrl(solver_approach,rec)
	% function return a list of ctrl structures and a string
	% to lable the approach.
	% solver_approach supported ::
	% 10 :: Based on primal Schur complement S=A+B1T*C^{-1} B2
	%       inverted with AGMG
	% 11 :: Based on dual Schur complement factorization with
	%       ~A=diag(A) ~S=-(C+B2*diag(A)^{-1} B1T)
	%       with the latter solved with AGMG
	% 13 :: HSS preconditioner

	
	ctrls=struct([]);
	approach_descriptions=[];

	if (solver_approach == 1)
		% approach based on direct solver
		ground=1;
		ground_node=1;
		reduced = 1;
		controls = struct('ground',ground,...
											'ground_node',ground_node,...
											'sol',solver_approach,...
											'reduced',reduced )
		approach_string='direct'

		ctrls=[ctrls,controls];
		approach_descriptions=[approach_descriptions,approach_string];
		

		
	elseif (solver_approach == 10)
		% approach based on inversion of Primal Schur Complement

		
		% set here bicgstab,gmres,fgmres (for non stationary prec)
		ctrl_outer = ctrl_solver;
		ctrl_outer.init('fgmres',1e-5,400,0.0,0);

		% preconditioner approach
		outer_prec='full';

		% reduce system to the primal schur complement
		solve_primal=0;


		% permute entries to minimize bandwidth S.
		% It will destroy its block structure.
		permute=0;

		
		% select assembly of S=A+Mtt+(Mtx+Mtx')+Mxx
		assembly_S_options=[...
												 "full",...
												 "A_Mtt",...
												 "harmonic",...
												 "A_Mtt_Mxx",...
												 "A_Mtt_Mxx_lamped",...
												 "A_Mtx_Mxt",...
												 "A_Mtt_Mtx_Mxt",...
												 "A_Mtt_Mtx_Mxt_blockdiagMxx"];
		for assembly_S = assembly_S_options([1])

			
			% cycle approach for S inversion
			% full:        :~S=S
			% block_triang :~S=upper block triangular of S
			% block_diag   :~S=block diagonal of S
			inverses11=["full","block_triang","block_diag"];
			for inverse11=inverses11(1)
				

				% scale
				diagonal_scaling=0;

				% handle singularity of S. 
				% ground==1 ground solution at zero in one node
				% ground==0 no grounding
				ground=0;
				ground_node=1;

				% S=S+relax*I			
				relax_inv11=1e-10;

				% controls for primal Schur complement
				ctrl_inner11=ctrl_solver;
				ctrl_inner11.init('agmg',1e-1,200,1.0,0,'agmg1e-1');

				% solver for C. It depends on the recostruction
				ctrl_inner22=ctrl_solver;
				if (rec==1)
					ctrl_inner22.init('diag',... %approach
														1e-6,... %tolerance
														100,...% itermax
														0.0,... %omega
														0,'diag'); %verbose
				else
					ctrl_inner22.init('agmg',... %approach
														1e-1,... %tolerance
														100,...% itermax
														0.0,... %omega
														0,'agmg1e-1'); %verbose
				end
				
				% store all controls
				controls = struct('ground',ground,...
													'ground_node',ground_node,...
													'sol',solver_approach,...
													'outer_prec',outer_prec,...
													'solve_primal',solve_primal,...
													'diagonal_scaling',diagonal_scaling,...
													'assembly_S',assembly_S,...
													'permute',permute,...
													'ctrl_outer',ctrl_outer,...
													'inverse11',inverse11,...
													'relax_inv11',relax_inv11,...
													'ctrl_inner11',ctrl_inner11,...
													'ctrl_inner22',ctrl_inner22);

				ctrls=[ctrls,controls];
				
				
				approach_string=strcat('SchurPrimal_',...
															 'ground=',num2str(ground),'_',...
															 'diagscal=',num2str(diagonal_scaling),'_',...
															 'outer=',ctrl_outer.approach,'_',...
															 'prec=',outer_prec,'_',...
															 'invS=',inverse11,'_',...
															 'solverS=',ctrl_inner11.label);

				approach_descriptions=[approach_descriptions,approach_string];
			end
		end
	elseif (solver_approach==11)
		% verbosity level in linear solver
		verbose = 0; 
		
		% set here fgmres (for non stationary prec), bicgstab,gmres, pcg
		outer_solvers={'bicgstab'  ,'gmres','fgmres' ,'pcg'};
		for isolver=[3]
			ctrl_outer=ctrl_solver;
			ctrl_outer.init(outer_solvers{isolver},1e-05,4000,0.0,0);

			% left or right preconditioner
			left_right='right';
			
			% external prec approach
			% full         :: P = Factorization in eq.
			% lower_triang :: P = ( A    )
			%                     ( B2 S )
			% lower_triang :: P = ( A  B1T )
			%                    (    S   )
			% identity     :: no preconditioner
			outer_precs= ["full" ,"lower_triang"  ,"upper_triang","identity"];
			for outer_prec=outer_precs([1,3])
				
				% relax A
				relax_inv11=0e-12;
				
				% set grounded_node>0 to gorund the potential in grounded node
				ground=0;
				diagonal_scaling=0;

				% set here list of solvers for block 22 
				solvers={'agmg'    ,'direct','krylov'  ,'incomplete'};
				iters  ={20,1      ,1        ,10       ,0           };
				label  ={'agmg1e-1','direct','krylov10','incomplete'};
				relax_inv22=0;				
				for isol=[1];%1:length(solvers)
					ctrl_inner22=ctrl_solver;
					ctrl_inner22.init(solvers{isol},1e-1,iters{isol},1.0,0,...
														sprintf('%s%1.1e',solvers{isol},1e-1));
					controls = struct('ground',ground,...
														'diagonal_scaling',diagonal_scaling,...
														'sol',solver_approach,...
														'outer_prec',outer_prec,...
														'left_right',left_right,...
														'ctrl_inner22',ctrl_inner22,...
														'ctrl_outer',ctrl_outer,...
														'verbose',verbose,...
														'relax_inv11',relax_inv11,...
														'relax_inv22',relax_inv22);

					ctrls=[ctrls,controls];
					

					approach_string=strcat('SIMPLE_','ground',num2str(ground),...
																 '_diagscal',num2str(diagonal_scaling),...
																 ctrl_outer.approach,'_',...
																 left_right,'_',outer_prec,'_prec_',...
																 'invSCA',ctrl_inner22.label);

					approach_descriptions=[approach_descriptions,approach_string];
				end
			end
		end
	elseif (solver_approach == 13)
		% global controls
		compute_eigen=0;
		verbose=0;

		
		ctrl_outer=ctrl_solver;
		% set here bicgstab,gmres,fgmres (for non stationary prec)
		ctrl_outer.init('fgmres',1e-5,3000);

		% left or right prec.
		left_right='right';

		% ground system
		ground=0;
		% diagonal scaling
		diagonal_scaling=1;


		% alpha relaxation
		alphas={0.5};
		for ll=1
			alpha=alphas{ll};

			ctrl_innerS=ctrl_solver;

			approach_inverse_A='block';
			%approach_inverse_A='full';

			% in S matrix invert w.r.t. the primal or the dual schur complement
			approaches_S={'primal','dual'};
			for ll=2
				approach_inverse_S=approaches_S{ll};

				% P=S*H or H*S
				approach_prec='SH';
				%approach_prec='HS';
				
				solversA={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete'};
				itersA  ={1      ,100       ,1       ,1        ,100        ,  1         };
				labelA  ={'agmg1','agmg1e-1' ,'direct','krylov1','krylov100','incomplete'};

				solversS={'agmg'  ,'agmg'  ,'direct','krylov' ,'krylov'  ,'incomplete','diag'};
				itersS  ={1      ,100       ,1       ,1        ,100        ,  1         ,1};
				labelS  ={'agmg1','agmg1e-1' ,'direct','krylov1','krylov100','incomplete','diag'};

				ctrl_innerC=ctrl_solver;
				%ctrl_innerC.init('diag',1e-1,100,1.0,0,'diag');
				ctrl_innerC.init('agmg',1e-1,100,1.0,0,'agmg1e-1');
				
				
				for isol=[2];%1:length(solversA)
					% set here other approximate inverse of A
					ctrl_innerA=ctrl_solver;
					ctrl_innerA.init(solversA{isol},1e-1,itersA{isol},1.0,0,labelA{isol});

					for j=[2];%1:length(solversS)
						% set here other approximate inverse of S
						
						ctrl_innerS.init(solversS{j},1e-1,itersS{j},1.0,0,labelS{j});


						
						controls = struct('ground',ground,...
															'sol',solver_approach,...
															'alpha',alpha,...
															'approach_prec',approach_prec,...
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

						approach_string=strcat('HSS_grounded',num2str(ground),...
																	 '_diagscal',num2str(diagonal_scaling),...
																	 approach_prec,'alpha',num2str(alpha),...
																	 '_invA',labelA{isol},...
																	 '_invC',ctrl_innerC.label,...
																	 '_invS_',approach_inverse_S,labelS{j});

						ctrls=[ctrls,controls];
						approach_descriptions=[approach_descriptions,approach_string];
						
					end
				end
			end
		end
	elseif (solver_approach == 20)
		% globals controls
		verbose=0;
		
		% list of outer solvers
		outer_solvers={'bicgstab'  ,'gmres','fgmres' ,'pcg'};
		
		% set here fgmres (for non stationary prec), bicgstab,gmres, pcg
		for isolver=[3]%1:length(outer_solvers)
			ctrl_outer=ctrl_solver;
			ctrl_outer.init(outer_solvers{isolver},1e-5,400,0.0,0);
			
			
			% external prec appraoch
			outer_precs={'full' ,'lower_triang'  ,'upper_triang','identity'};
			nouter_precs=length(outer_precs);

			% left or right preconditioner
			left_right='left';
			% fgmres works only with right prec
			if (strcmp(outer_solvers,'fgmres'))
				left_right='right';
			end

			% node for grounding
			grounded_node=0;
			
			diagonal_scaling=0;
			
			
			for iprec=[3]%1:nouter_precs
				outer_prec=outer_precs{iprec};
				
				% inverse approach
				%inverse11='full';
				inverse11='diag';
				
				% relaxation
				% A = A + relax * Id
				relax4inv11=1e-12;

				for precision_primal=1
					% set here other approximate inverse of block11
					ctrl_inner11=ctrl_solver;
					ctrl_inner11.init('agmg',... %approach
														10^(-precision_primal),... %tolerance
														20,...% itermax
														0.0,... %omega
														0,strcat('agmg',num2str(precision_primal))); %verbose
					
					% possible choices for dual schur complement
					inverses22={'diagA','perm','lsc','full','commute','commute_bis','commute_ter','commute_quater'};
					null_space=0;
					lrb='r';
					ground=0;

					% for commute approach
					mode_inverse22=1;
					
					
					for ia = [8];
						inverse22=inverses22{ia};
						if ( contains(inverse22,'commute') )
							label_inverse22=strcat(inverse22,'null_space',num2str(null_space),'lrb', lrb);
						else
							label_inverse22=strcat(inverse22);
						end
						
						% set here list of solvers for block 22 
						solvers={'agmg' ,'agmg'  ,'agmg' ,'direct','krylov' ,'krylov'  ,'incomplete','diag','krylov_no_prec'};
						iters  ={1      ,10      ,400,1       ,1        ,10        ,  1          ,0,10};
						label  ={'agmg1','agmg10','agmg1e-1','direct','krylov1','krylov10','incomplete','diag','purekrylov'};
						relax4inv22=1e-12;

						tolerance_preprocess=1e-2;

						for precision_dual=1
							for isolver=[3];%1:length(solvers)
								ctrl_inner22=ctrl_solver;
								ctrl_inner22.init('agmg',10^(-precision_dual),30,0.0,0,strcat('agmg',num2str(precision_dual)));%,...
								%sprintf('%s%1.1e',solvers{i},dual_schur_tol));
								controls = struct('indc',grounded_node,...
																	'ground',ground,...
																	'sol',solver_approach,...
																	'null_space',null_space,...
																	'lrb',lrb,...
																	'diagonal_scaling',diagonal_scaling,...
																	'outer_prec',outer_prec,...
																	'left_right',left_right,...
																	'ctrl_inner11',ctrl_inner11,...
																	'ctrl_inner22',ctrl_inner22,...
																	'ctrl_outer',ctrl_outer,...
																	'verbose',verbose,...
																	'inverse11',inverse11,...
																	'relax4inv11',relax4inv11,...
																	'inverse22',inverse22,...
																	'relax4inv22',relax4inv22,...
																	'mode_inverse22',mode_inverse22,...
																	'tolerance_preprocess',tolerance_preprocess);
								
								
								approach_string=strcat('augmented_',...
																			 'ground',num2str(ground),'_',...
																			 'outer_prec',outer_prec,'_',...
																			 'invA',ctrl_inner11.label,'_',...
																			 'inverse22',label_inverse22,'_',...
																			 'mode',num2str(mode_inverse22),'_',...
																			 'invSCA',ctrl_inner22.label);
								
								ctrls=[ctrls,controls];
								approach_descriptions=[approach_descriptions,approach_string];
							end
						end
					end
				end
			end
		end
	elseif (solver_approach==220)
		% global controls
		verbose=0;
		
		% set here fgmres (for non stationary prec), bicgstab,gmres, pcg
		ctrl_outer=ctrl_solver;
		ctrl_outer.init('fgmres',1e-05,400,0.0,0);
		left_right='right';
		
		% external prec appraoch
		outer_precs=["full","lower_triang","upper_triang"];
		for outer_prec=outer_precs([1])
			% set grounded_node>0 to gorund the potential in grounded node
			grounded_node=0;
			diagonal_scaling=0;

			% primal is based on inverse of J=  [ A B1T ] and inverse of Dr
			%                                   [ B2 -C ]
			% dual  is based on inverse of  J = [ A B1T ] and inverse of Dr+Ds S^{-1} M
			%                                   [ B2 0  ]
			%approach_Schur_slack='dual';
			approach_Schur_slack='primal';


			% solving block J
			outer_prec_J='upper_triang'; % we invert A just once
			
			% inverse A with its diagonal
			relax4_inv11=1e-12;
			ctrl_inner11=ctrl_solver;
			for precions_primal=[5e-1,1e-1];
				ctrl_inner11.init('agmg',precions_primal,10,1.0,0);

				approaches_Schur_rho=[...
															 "commute_Btilde",...
															 "commute_right",...
															 "commute_left",...
															 "commute_both",...
															 "lsc",...
															 "diag"];
				for approach_Schur_rho = approaches_Schur_rho(1)
					% inverse schur_rho 
					relax4_inv22=1e-12;
					ctrl_inner22=ctrl_solver;
					ctrl_inner22.init('agmg',1e-1,20,1.0,0);

					
					
					% inverse schur_s 
					relax4_inv33=1e-12;
					ctrl_inner33=ctrl_solver;
					ctrl_inner33.init('diag',1e-1,10,1.0,0);
					%ctrl_inner33.init('bicgstab',1e-1,100,1.0,1,'S33');


					controls = struct('indc',grounded_node,...
														'diagonal_scaling',diagonal_scaling,...
														'sol',solver_approach,...
														'outer_prec',outer_prec,...
														'outer_prec_J',outer_prec_J,...
														'left_right',left_right,...
														'approach_Schur_rho',approach_Schur_rho,...
														'relax4inv11',relax4_inv11,...
														'ctrl_inner11',ctrl_inner11,...
														'relax4inv22',relax4_inv22,...
														'ctrl_inner22',ctrl_inner22,...
														'approach_Schur_slack',approach_Schur_slack,...
														'relax4inv33',relax4_inv33,...
														'ctrl_inner33',ctrl_inner33,...											
														'ctrl_outer',ctrl_outer,...
														'verbose',verbose);
					

					approach_string=strcat('commute_ground',num2str(grounded_node),'_',...
																 'diagscal',num2str(diagonal_scaling),'_',...
																 'outer=',approach_Schur_slack,'_',...
																 'outprec=',outer_prec,'_',...
																 'outer_prec_J=',outer_prec_J,'_',...
																 'invS=',approach_Schur_rho,'_',...
																 'invA',ctrl_inner11.label,'_',...
																 'inv22',ctrl_inner22.label,'_',...
																 'inv33',ctrl_inner33.label);

					ctrls=[ctrls,controls];
					approach_descriptions=[approach_descriptions,approach_string];
					
					
				end
			end
		end
	end
end
