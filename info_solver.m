classdef info_solver <handle
   properties
     flag;
     res;
     realres;
     iter;
     rhsnorm;
     approach_used;
   end
   methods
     function obj = print(obj)
       fprintf('%s flag=%d iter=%d realres=%9.4e res=%9.4e rhs=%9.4e\n',...
	       obj.approach_used,...
	       obj.flag,...
	       obj.iter,...
	       obj.realres,...
	       obj.res,...
	       obj.rhsnorm)
     end
   end
end
