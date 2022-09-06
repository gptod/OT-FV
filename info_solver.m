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
     sumsol;
   end
   methods
     function obj = print(obj,fID);
       if (~exist('fid','var') )
	 fid=1;
       end
       fprintf(fid,'%s %s flag=%d it=%3d rel.res=%1.1e |res|=%1.1e sum(sol)=%1.1e |rhs|=%1.1e sum(rhs)=%1.1e\n',...
	       obj.label,...
	       obj.approach_used,...
	       obj.flag,...
	       uint64(obj.iter),...
	       obj.realres,...
	       obj.res,obj.sumsol,...
	       obj.rhsnorm,obj.balance);
     end
   end
end
