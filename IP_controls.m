% class spatial information of grid
classdef IP_controls <handle
  properties
		eps_0  = 1e-5; % tolerance for the IPM
		eps_mu = 1e-6; % tolerance Newton
		k2max = 15; % maximum number of inner (Newton) iterations
		k1max = 13; % maximum number of outer iterations
		theta0 = 0.2; % decay ratio for the perturbation parameter mu
		theta_min = 0.2;
		theta_max = 0.2;
		alpha_min = 0.1; % minimal step length accepted
		mu0_up = 0;
		mu=1.0; %relaxation parameters

		%% mode for linear soslver setting
		% kel = 0 fixed tolerance
		% kel = 1 Kelly on
		kel  = 0; 
		
		min_outer_tol= 1e-5; %lower bound for solving linear system
		
		% verbosity
		verbose=2;
		
		% save hystory for save into file_log hystory of algorithm
		save_log=0; % on/off control
		file_log='l2opt_IP.log' ; % full path destionation
		% save file_csv wiht info on IP and linear algebra solvers
		save_csv=0; % on/off control
		file_csv='l2opt_IP.csv'; % full path destionation
		% save file into h5 format
		save_h5=0; % 0 : do not save; 1: save results results IP cycle; 2: save inner iterations
		file_h5='l2opt_IP.h5'; % full path directory
	end
  methods
		% constructor
    function obj = init()
    end
   	% destructor
		function obj = kill(obj)
		end

		% info
		function obj = info(obj,fid)
      if (~exist('fid','var') )
				fid=1;
			end
		end
	end
end
