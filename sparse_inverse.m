classdef sparse_inverse <handle
   properties
     matrix;
     nequ;
     is_symmetric;
     ctrl;
     info;
     % agmg solver
     mgsolver;
     % incomplete factor in LU factorization A~IL*IU
     IL; 
     IU;
     % M, N matrix for Gauss-Seidel and Jacobi
     M; 
     N; 
   end
   methods
     function obj = init(obj,matrix,ctrl)
       obj.matrix=matrix;
       obj.nequ=size(matrix,1);
       obj.is_symmetric = norm(nonzeros(matrix)-nonzeros(matrix'))<1e-14;
       obj.ctrl=ctrl;
       obj.info=info_solver;

       if ( strcmp(ctrl.approach,'\'))
       elseif ( strcmp(ctrl.approach,'agmg') )
	  if (obj.is_symmetric)
	    icg=1;
	 else
	   icg=obj.ctrl.nrestart;
	 end
	 %obj.mgsolver=agmg(matrix,[],icg,[],ctrl.itermax,0,[],1);
       elseif ( strcmp(ctrl.approach,'krylov'))
	 if ( obj.is_symmetric)
	   obj.IL = ichol(matrix, struct('type','ict','droptol',1e-3));
	 else
	   setup.type = 'crout';
	   setup.milu = 'row';
	   setup.droptol = 1e-5;
	   [obj.IL,obj.IU] = ilu(matrix,setup);
	 end
       elseif( strcmp(ctrl.approach , 'incomplete') )
	 if ( obj.is_symmetric)
	   obj.L = ichol(matrix, struct('type','ict','droptol',1e-3));
	 else
	   setup.type = 'crout';
	   setup.milu = 'row';
	   setup.droptol = 1e-5;
	   [obj.L,obj.U] = ilu(matrix,setup);
	 end
       elseif ( strcmp(obj.ctrl.approach, 'gs') )
	 obj.D=spdiags(diag(obj.matrix),0,obj.nequ,obj.nequ);
	 obj.M=(1/omega)*obj.D+tril(obj.matrix,-1);
	 obj.N=((1-omega)/omega)*obj.D-triu(obj.matrix,1);
       elseif ( strcmp(obj.ctrl.approach, 'jacobi') )
	 obj.D=spdiags(diag(obj.matrix),0,obj.nequ,obj.nequ);
	 obj.M=obj.D;
	 obj.N=obj.M-obj.matrix;
       end
       
     end
     % sol ~= A^{-1}(rhs) 
     function sol = apply(obj,rhs,initial_guess)
       if (~exist('initial_guess','var'))
	 initial_guess=zeros(obj.nequ,1);
       end
		 
       if ( strcmp(obj.ctrl.approach ,'\'))
	 sol=obj.matrix\rhs;
	 obj.info.approach_used = 'backslash';
       elseif ( strcmp(obj.ctrl.approach ,'agmg'))
	 if (obj.is_symmetric)
	   icg=1;
	 else
	   icg=obj.ctrl.nrestart;
	 end
	 [sol,obj.info.flag, obj.info.res, obj.info.iter,obj.info.resvec]=...
	 agmg(obj.matrix,rhs,icg,obj.ctrl.tolerance,obj.ctrl.itermax,0,initial_guess,0);	 
	 obj.info.approach_used = 'agmg';
       elseif( strcmp(obj.ctrl.approach, 'krylov'))
	 if ( obj.is_symmetric )
	   [sol,obj.info.flag,obj.info.res,obj.info.iter,obj.info.resvec]=...
           pcg(obj.matrix,rhs,...
               obj.ctrl.tolerance,obj.ctrl.itermax,[],@(x) obj.IL\(obj.IL'\x),[],initial_guess);
	   obj.info.approach_used = 'pcg+ic';
	 else
	   [sol,obj.info.flag,obj.info.res,obj.info.iter,obj.info.resvec]=...
           bicgstab(obj.matrix,rhs,...
		 obj.ctrl.tolerance,obj.ctrl.itermax,[],@(x) obj.IU\(obj.IL\x),[],initial_guess);
	   obj.info.approach_used = 'bicgtab+ilu';
	 end

	%
	%  invert with incomplete factorization
	% 
       elseif( strcmp(obj.ctrl.approach , 'incomplete'))
	 if ( obj.is_symmetric )
	   sol=obj.IL\(obj.IL'\rhs);
	   obj.info.approach_used = 'ic';
	 else
	   sol=obj.IL\(obj.IU\rhs);
	   obj.info.approach_used = 'ilu';
	 end
       elseif ( strcmp(obj.ctrl.approach, 'gs') || strcmp(obj.ctrl.approach, 'jacobi')  )
	   matrix version
	   obj.info.iter=0;
	   sol=initial_guess;
	   res=norm(obj.matrix*sol-rhs);
	   resv=res;
	   omega=ctrl_solver.omega
	   obj.info.iter=0
	   while res > obj.ctrl.tolerance && k<obj.ctrl.itermax
	     obj.info.iter=obj.info.iter+1;
	     sol=obj.M\(obj.N*sol+rhs);
	     res=norm(obj.matri*x-rhs);
	     resv=[resv;res];
	   end
	   obj.info.approach_used = 'gs';
	   obj.info.resvec = resvec;
       end
       obj.info.rhsnorm=norm(rhs);
       obj.info.realres=norm(obj.matrix*sol-rhs)/norm(rhs);
     end
     % destructor
     function obj = kill(obj)
       if ( strcmp(obj.ctrl.approach ,'agmg'))
	 %obj.mgsolver=agmg(obj.matrix,[],1,[],1000,0,[],-1);
       end
       clear obj.matrix;
       clear obj.is_symmetric;
       clear obj.ctrl;
       clear obj.info;
       clear obj.IL;
       clear obj.IU;
     end
   end
end
