function [] = plot_rhos(grid_rho,rho_in,rho,rho_f)
	ncell_rho = grid_rho.ncell;
	N = size(rho,1) / ncell_rho;
	figure
	fig=patch('Faces',grid_rho.cells(:,2:end),...
						'Vertices',grid_rho.nodes,...
						'edgecolor','none',...
						'FaceVertexCData',rho_in,'FaceColor','flat');
	colorbar
	axis square
	axis off
	caxis([0 max(rho)])
	%str =num2str(0);
	%outname=strcat('cross_harm_rhok',str,'.jpg');
	%saveas(fig,outname)
	figure
	fig=patch('Faces',grid_rho.cells(:,2:end),'Vertices',grid_rho.nodes,...
						'edgecolor','none','FaceVertexCData',rho_f,'FaceColor','flat');
	colorbar
	axis square
	axis off
	caxis([0 max(rho)])
	%str =num2str(nk+1);
	%outname=strcat('cross_harm_rhok',str,'.jpg');
	%saveas(fig,outname)
	%dsp = 2;
	dsp = ceil(N/2);
	for k=1:N
		%k=4;
		rhok = rho((k-1)*ncell_rho+1:k*ncell_rho);
		figure
		fig=patch('Faces',grid_rho.cells(:,2:end),...
							'Vertices',grid_rho.nodes,...
							'edgecolor','none','FaceVertexCData',rhok,'FaceColor','flat');
		colorbar
		axis square
		axis off
		caxis([0 max(rho)])
		str =num2str(k);
		outname=strcat('gauss_linear_rhok',str,'.jpg');
		saveas(fig,outname)
	end

end
