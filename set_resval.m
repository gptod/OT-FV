classdef set_resval <handle
    properties
        theta=0.9
        etamax=0.1
        etak
        fnrm
        fnrm0
        stop_tol
    end
    methods
        function obj=set(obj,kel,eps,varargin)
            if kel==0
                obj.etak = eps;
            else
                if isempty(obj.fnrm0)
                    delta = varargin{1};
                    obj.theta=.9;
                    %etamax=.9;
                    obj.etamax = 0.1;
                    obj.etak = obj.etamax;
                    obj.fnrm = delta;%/sqrt(tnp+2*tnr2h);
                    obj.fnrm0 = 1;
                    %atol = 1e-6; rtol = 1e-8;
                    %atol = eps; rtol = 0.01*eps;
                    %obj.stop_tol = atol + rtol*obj.fnrm;
                    obj.stop_tol = eps;
                else
                    delta = varargin{1};
                    obj.fnrm = delta;
                    rat=obj.fnrm/obj.fnrm0;
                    obj.fnrm0=obj.fnrm;
                    % adjust eta:
                    if obj.etamax > 0
                        etaold = obj.etak;
                        etanew = obj.theta*rat*rat;
                        if obj.theta*etaold*etaold > .1
                            etanew = max(etanew,obj.theta*etaold*etaold);
                        else
                            obj.etak = min([etanew,obj.etamax]);
                        end
                        %obj.etak = min(obj.etamax,max(obj.etak,0.5*obj.stop_tol/obj.fnrm));
                        obj.etak = min(obj.etamax,max(obj.etak,0.5*obj.stop_tol));
                    end
                end
            end
        end
    end
end

