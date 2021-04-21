classdef info_solver <handle
   properties
     flag;
     res;
     realres;
     resvec;
     iter;
     rhsnorm;
     approach_used;
   end
   methods
     function obj = print(obj,fID);
       if (~exist('fid','var') )
	 fid=1;
       end
       fprintf(fid,'%s flag=%d iter=%3d res=%9.4e rhs=%9.4e\n',...
	       obj.approach_used,...
	       obj.flag,...
	       uint8(obj.iter),...
	       obj.res,...
	       obj.rhsnorm);
     end
   end
end
