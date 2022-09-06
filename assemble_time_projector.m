function Pt = assemble_time_projector(N,ncell)
	
	Pt = 1 * (spdiags(ones((N)*ncell,1),0,(N)*ncell,(N+1)*ncell)) + ...
			 0 * (spdiags(ones((N)*ncell,1),ncell,(N)*ncell,(N+1)*ncell));

end
