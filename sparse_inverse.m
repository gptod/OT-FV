classdef sparse_inverse <handle
   properties
     matrix;
     is_symmetric;
     ctrl;
     info;
     mgsolver;
     IL;
     IU;
   end
   methods
     function obj = init(obj,matrix,ctrl)
       obj.matrix=matrix;
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
       end
       
     end
     % sol ~= A^{-1}(rhs) 
     function sol = apply(obj,rhs)
       if ( strcmp(obj.ctrl.approach ,'\'))
	 sol=obj.matrix\rhs;
	 obj.info.approach_used = 'backslash';
       elseif ( strcmp(obj.ctrl.approach ,'agmg'))
	 if (obj.is_symmetric)
	   icg=1;
	 else
	   icg=obj.ctrl.nrestart;
	 end
	 [sol,obj.info.flag, obj.info.res, obj.info.iter]=...
	 agmg(obj.matrix,rhs,icg,obj.ctrl.tolerance,obj.ctrl.itermax,0,[],0);	 
	 obj.info.approach_used = 'AGMG';
       elseif( strcmp(obj.ctrl.approach, 'krylov'))
	 if ( obj.is_symmetric )
	   [sol,obj.info.flag,obj.info.res,obj.info.iter]=...
           pcg(obj.matrix,rhs,...
               obj.ctrl.tolerance,obj.ctrl.itermax,[],@(x) obj.IL\(obj.IL'\x));
	   obj.info.approach_used = 'PCG+IC';
	 else
	   [sol,obj.info.flag,obj.info.res,obj.info.iter]=...
           bicgstab(obj.matrix,rhs,...
		 obj.ctrl.tolerance,obj.ctrl.itermax,[],@(x) obj.IU\(obj.IL\x));
	   obj.info.approach_used = 'BICGTAB+ILU';
	 end
	 
       elseif( strcmp(obj.ctrl.approach , 'incomplete'))
	 if ( obj.is_symmetric )
	   sol=obj.IL\(obj.IL'\rhs);
	   obj.info.approach_used = 'IC';
	 else
	   sol=obj.IL\(obj.IU\rhs);
	   obj.info.approach_used = 'ILU';
	 end
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
