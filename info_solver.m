classdef info_solver <handle
  properties
    label='';
     flag;
     res;
     realres;
     resvec;
     iter;
     rhsnorm;
     balance;
     approach_used;
   end
   methods
     function obj = print(obj,fID);
       if (~exist('fid','var') )
	 fid=1;
       end
       fprintf(fid,'%s %s flag=%d iter=%3d res=%9.4e realres=%9.4e rhs=%9.4e sum(rhs)=%9.4e\n',...
	       obj.label,...
	       obj.approach_used,...
	       obj.flag,...
	       uint64(obj.iter),...
	       obj.res,...
	       obj.realres,...
	       obj.rhsnorm,obj.balance);
     end
   end
end
