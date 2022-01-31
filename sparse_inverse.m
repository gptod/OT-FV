classdef sparse_inverse <handle
  properties
    dimblock=0;
    name='empty';
    preprocess_agmg=0;
    matrix;
    inverse_matrix_diagonal;
    nequ;
    is_symmetric;
    ctrl;
    info_inverse;
    % dense factorizations
    matrix_decomposed; 
    % agmg solver
    mgsolver;
    % incomplete factor in LU factorization A~IL*IU
    IL; 
    IU;
    % M, N matrix for Gauss-Seidel and Jacobi
    M; 
    N;
    cumulative_iter=0;
    cumulative_application=0;
    init_cpu=0;
    cumulative_cpu=0; 
   end
   methods
     function obj = init(obj,matrix,ctrl)
       obj.init_cpu=0;
       init_cpu=tic;
       obj.matrix=matrix;
       obj.nequ=size(matrix,1);
       obj.is_symmetric = norm(nonzeros(matrix)-nonzeros(matrix'))<1e-12;
       obj.ctrl=ctrl;
       obj.info_inverse=info_solver;
       obj.cumulative_iter=0;
       obj.cumulative_cpu=0;
       obj.cumulative_application=0;
       obj.preprocess_agmg=ctrl.preprocess_agmg;

       if ( contains(ctrl.approach,'direct'))
	 obj.matrix_decomposed = decomposition(obj.matrix); 
       elseif ( contains(ctrl.approach,'agmg') )
	  if (obj.is_symmetric)
	    icg=1;
	 else
	   icg=obj.ctrl.nrestart;
	  end
	  if ( obj.preprocess_agmg > 0 )
	    agmg_str=sprintf('agmg%d',obj.preprocess_agmg);
				%disp( strcat('init',agmg_str));
	    verbose_agmg=0;
	    obj.mgsolver=feval( agmg_str, obj.matrix,[],icg,verbose_agmg,ctrl.itermax,0,[],1);
	  end
	  
       elseif ( strcmp(ctrl.approach,'diag') )
	 obj.inverse_matrix_diagonal= 1.0./spdiags(matrix,0);
       elseif ( strcmp(ctrl.approach,'krylov'))
	 if ( obj.is_symmetric)
	   obj.IL = ichol(matrix, struct('type','ict','droptol',1e-3));
	 else
	   setup.type = 'crout';
	   setup.milu = 'row';
	   setup.droptol = 1e-3;
	   [obj.IL,obj.IU] = ilu(matrix,setup);
	 end
       elseif ( strcmp(ctrl.approach,'krylov_no_prec'))
	 % nothing to do
       elseif( strcmp(obj.ctrl.approach , 'incomplete') )
	 if ( obj.is_symmetric)
	   obj.IL = ichol(matrix, struct('type','ict','droptol',1e-5));
	 else
	   setup.type = 'crout';
	   setup.milu = 'row';
	   setup.droptol = 1e-5;
	   [obj.IL,obj.IU] = ilu(matrix,setup);
	 end
       elseif ( strcmp(obj.ctrl.approach, 'gs') )
	 obj.D=spdiags(diag(obj.matrix),0,obj.nequ,obj.nequ);
	 obj.M=(1/omega)*obj.D+tril(obj.matrix,-1);
	 obj.N=((1-omega)/omega)*obj.D-triu(obj.matrix,1);
       elseif ( strcmp(obj.ctrl.approach, 'jacobi') )
	 obj.D=spdiags(diag(obj.matrix),0,obj.nequ,obj.nequ);
	 obj.M=obj.D;
	 obj.N=obj.M-obj.matrix;
       else
	 fprintf('In sparse_inverse.init: approach %s not supported. Execution will be stopped\n',obj.ctrl.approach)
	 return
       end
       obj.init_cpu=toc(init_cpu);
       
     end
     % sol ~= A^{-1}(rhs) 
     function sol = apply(obj,rhs,initial_guess)
       if (~exist('initial_guess','var'))
	 initial_guess=zeros(obj.nequ,1);
       end

       result = sum(isnan(rhs(:)));
       if (result>0)
	 disp('Nan in RHS')
       end
       
       apply_cpu=tic;
       if ( strcmp(obj.ctrl.approach ,'direct'))
	 sol= obj.matrix_decomposed\rhs;
	 obj.info_inverse.approach_used = 'direct';
	 obj.info_inverse.iter=0;
	 obj.info_inverse.flag=0;
       elseif ( contains(obj.ctrl.approach ,'agmg'))
	 if (obj.is_symmetric)
	   icg=1;
	 else
	   icg=obj.ctrl.nrestart;
	 end
	 verbose_agmg=0;
	 if (obj.preprocess_agmg>0)
	   if (strcmp(obj.ctrl.approach ,'agmg'))
	     jobagmg=2;
	     agmg_str=sprintf('agmg%d',obj.preprocess_agmg);
	     [sol,obj.info_inverse.flag, obj.info_inverse.res, obj.info_inverse.iter,obj.info_inverse.resvec]=...
	     feval( agmg_str, obj.matrix,rhs,icg,obj.ctrl.tolerance,obj.ctrl.itermax,verbose_agmg,[],jobagmg);
	     obj.info_inverse.approach_used = 'agmg';
	   elseif(strcmp(obj.ctrl.approach ,'precagmg'))
	     jobagmg=3;
	     agmg_str=sprintf('agmg%d',obj.preprocess_agmg);
	     [sol,obj.info_inverse.flag, obj.info_inverse.res, obj.info_inverse.iter,obj.info_inverse.resvec]=...
	     feval( agmg_str, obj.matrix,rhs,icg,[],[],verbose_agmg,[],jobagmg);
	     obj.info_inverse.approach_used = 'precagmg';
	   end
	 else
	   jobagmg=0;
	   [sol,obj.info_inverse.flag, obj.info_inverse.res, obj.info_inverse.iter,obj.info_inverse.resvec]=...
	   agmg(obj.matrix,rhs,icg,obj.ctrl.tolerance,obj.ctrl.itermax,verbose_agmg,initial_guess,jobagmg);
	   obj.info_inverse.approach_used = 'agmg';
	 end
	 
       elseif ( strcmp(obj.ctrl.approach,'diag') )
	 sol=obj.inverse_matrix_diagonal.*rhs;
	 obj.info_inverse.iter=1;
	 obj.info_inverse.flag=0; 
       elseif( strcmp(obj.ctrl.approach, 'krylov'))
	 if ( obj.is_symmetric )
	   [sol,obj.info_inverse.flag,obj.info_inverse.res,obj.info_inverse.iter,obj.info_inverse.resvec]=...
           pcg(obj.matrix,rhs,...
               obj.ctrl.tolerance,obj.ctrl.itermax, obj.IL, obj.IL',initial_guess);
	   obj.info_inverse.approach_used = 'pcg+ic';
	 else
	   [sol,obj.info_inverse.flag,obj.info_inverse.res,obj.info_inverse.iter,obj.info_inverse.resvec]=...
           bicgstab(obj.matrix,rhs,...
		 obj.ctrl.tolerance,obj.ctrl.itermax,obj.IL,obj.IU,initial_guess);
	   obj.info_inverse.approach_used = 'bicgtab+ilu';
	 end
       elseif( strcmp(obj.ctrl.approach, 'krylov_no_prec'))
	 if ( obj.is_symmetric )
	   [sol,obj.info_inverse.flag,obj.info_inverse.res,obj.info_inverse.iter,obj.info_inverse.resvec]=...
           pcg(obj.matrix,rhs,...
               obj.ctrl.tolerance,obj.ctrl.itermax, [], [],initial_guess);
	   obj.info_inverse.approach_used = 'pcg';
	 else
	   [sol,obj.info_inverse.flag,obj.info_inverse.res,obj.info_inverse.iter,obj.info_inverse.resvec]=...
           bicgstab(obj.matrix,rhs,...
		 obj.ctrl.tolerance,obj.ctrl.itermax,[],[],initial_guess);
	   obj.info_inverse.approach_used = 'bicgtab';
	 end

	%
	%  invert with incomplete factorization
	% 
       elseif( strcmp(obj.ctrl.approach , 'incomplete'))
	 if ( obj.is_symmetric )
	   sol=obj.IL\(obj.IL'\rhs);
	   obj.info_inverse.approach_used = 'ic';
	 else
	   sol=obj.IL\(obj.IU\rhs);
	   obj.info_inverse.approach_used = 'ilu';
	 end
	 obj.info_inverse.iter=0;
       elseif ( strcmp(obj.ctrl.approach, 'gs') || strcmp(obj.ctrl.approach, 'jacobi')  )
	   matrix version
	   obj.info_inverse.iter=0;
	   sol=initial_guess;obj.info_inverse.iter=0;
	   res=norm(obj.matrix*sol-rhs);
	   resv=res;
	   omega=ctrl_solver.omega
	   obj.info_inverse.iter=0
	   while res > obj.ctrl.tolerance && k<obj.ctrl.itermax
	     obj.info_inverse.iter=obj.info_inverse.iter+1;
	     sol=obj.M\(obj.N*sol+rhs);
	     res=norm(obj.matri*x-rhs);
	     resv=[resv;res];
	   end
	   obj.info_inverse.approach_used = 'gs';
	   obj.info_inverse.resvec = resvec;
       end
       cpu=toc(apply_cpu);

       obj.cumulative_application=obj.cumulative_application +1;
       obj.info_inverse.rhsnorm=norm(rhs);
       obj.info_inverse.balance=sum(rhs);
       obj.info_inverse.res=norm(obj.matrix*sol-rhs);
			 
       obj.info_inverse.realres=obj.info_inverse.res/norm(rhs);
       obj.cumulative_iter=obj.cumulative_iter+obj.info_inverse.iter;
       obj.cumulative_cpu=obj.cumulative_cpu+cpu;

       obj.info_inverse.label=obj.ctrl.label;
       obj.info_inverse.sumsol=sum(sol);
       if (obj.ctrl.verbose)
				 obj.info_inverse.print(obj.ctrl.unit_out);
       end
       
       
       if ( obj.info_inverse.flag ~= 0 )
	 obj.info_inverse();
	 %error('INNER ERROR FOR LINEAR SOLVER');
       end
     end
     % destructor
     function obj = kill(obj)
       if ( contains(obj.ctrl.approach ,'agmg'))
	 if( obj.preprocess_agmg>0 )
	   
	   agmg_str=sprintf('agmg%d',obj.preprocess_agmg);
	   %disp( strcat('killing',agmg_str));
	   z=...
	   feval( agmg_str, obj.matrix,[],1,[],1000,0,[],-1);
	 end
       end
       clear obj.matrix_decomposed;
       clear obj.matrix;
       clear obj.is_symmetric;
       clear obj.ctrl;
       clear obj.info_inverse;
       clear obj.IL;
       clear obj.IU;
     end

     % info
     function obj = info(obj,fid)
        if (~exist('fid','var') )
	  fid=1;
       end
       obj.ctrl.info(fid);
       if ( strcmp(obj.ctrl.approach ,'krylov'))
	 if ( obj.is_symmetric) 
	   fprintf(fid,'nnz(IL)/nnz(A),nnz(A) %f3 %d \n',...
		   nnz(obj.IL)/nnz(obj.matrix),nnz(obj.matrix));
	 else
	   fprintf(fid,'(nnz(IL),nnz(IU))/nnz(A)| nnz(A) %f3 %f3 %d\n',...
		   nnz(obj.IL)/nnz(obj.matrix),...
		   nnz(obj.IU)/nnz(obj.matrix),...
		   nnz(obj.matrix));
	 end
       end
     end
   end
end
