classdef ctrl_solver <handle
   properties
     approach;
     tolerance;
     itermax;
     nrestart=20;
     omega=1.0;
     verbose=0;
     preprocess_agmg=0;
     label='';
     unit_out=1;
   end
   methods
     function obj = init(obj,approach,tolerance,itermax,omega,verbose,label,preprocess_agmg,unit_out)
       obj.approach=approach;
       obj.tolerance=tolerance;
       obj.itermax=itermax;
       
       if (~exist('omega','var') )
				 obj.omega=1;
       else
				 obj.omega=omega;
       end
       if (~exist('verbose','var') )
				 obj.verbose=0;
       else
				 obj.verbose=verbose;
       end
       if (~exist('label','var') )
				 obj.label=sprintf('%s%1.1e',approach,tolerance);
				 if ( strcmp(obj.approach,'diag'))
					 obj.label='diag';
				 end
				 if strcmp( obj.approach , 'direct')
					 obj.label='direct';
				 end
			 else
				 obj.label=label;
			 end
		 

       if (~exist('preprocess_agmg','var') )
	 obj.preprocess_agmg=0;
       else
	 obj.preprocess_agmg=preprocess_agmg;
       end

       if (~exist('unit_out','var') )
	 obj.unit_out=1;
       else
	 obj.unit_out=unit_out;
       end
       
     end
     function obj = info(obj,fid)
        if (~exist('fid','var') )
	  fid=1;
       end
       fprintf(fid,'approach=%s tolerance=%9.3e itermax=%d \n',...
	       obj.approach,obj.tolerance,obj.itermax);
     end
   end
end

%ctrl_test=ctrl_solver;
%ctrl_test.init('agmg',1e-5,200);
%cltr_test.info()

% classdef info_solver
%    properties
%       flag
%       res
%       iter
%    end
% end
      

% classdef spinverse
%    properties
%       matrix
%       ctrl
%       info
%    end
%    methods
%       function obj = set.inverse(obj)
%          if (strcmpi(material,'aluminum') ||...
%                strcmpi(material,'stainless steel') ||...
%                strcmpi(material,'carbon steel'))
%             obj.Material = material;
%          else
%             error('Invalid Material')
%          end
%       end
%    end
% end

% function [sol] = define_inverse(matrix, rhs,ctrl)
%   if ( ctrl.approach == '\')
%     sol=matrix\rhs;
%   else if( ctrl.approach == 'agmg')
% 	 agmg(S,rhs,0,1e-6,200,0,[],0)
