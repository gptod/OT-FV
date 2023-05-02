function [] = plot_rhos(grid_rho,rho_in,rho,rho_f,dsp)
	dim =  size(grid_rho.nodes,2);
	if dim == 2
		ncells_rho = grid_rho.ncells;
		N = size(rho,1) / ncells_rho;
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
		%dsp = ceil(N/2);
		for k=dsp
			rhok = rho((k-1)*ncells_rho+1:k*ncells_rho);
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
	else
		ncells_rho = grid_rho.ncells;

    figure	
    Frho = scatteredInterpolant(grid_rho.cc,rho_in);
    Frho.ExtrapolationMethod = 'nearest';
    [X,Y,Z]=meshgrid(linspace(0,1,ceil((grid_rho.ncells)^(1/3))));
    V = Frho(X,Y,Z);
    xslice = 0.3;   
    yslice = [];
    zslice = 0.3;
    slice(X,Y,Z,V,xslice,yslice,zslice,'nearest')
    axis square
    colorbar
    
    figure	
    Frho = scatteredInterpolant(grid_rho.cc,rho_f);
    Frho.ExtrapolationMethod = 'nearest';
    [X,Y,Z]=meshgrid(linspace(0,1,ceil((grid_rho.ncells)^(1/3))));
    V = Frho(X,Y,Z);
    xslice = 0.7;   
    yslice = [];
    zslice = 0.7;
    slice(X,Y,Z,V,xslice,yslice,zslice,'nearest')
    axis square
    colorbar
    
		for k=dsp
			rhok = rho((k-1)*ncells_rho+1:k*ncells_rho);
			figure	
      Frho = scatteredInterpolant(grid_rho.cc,rhok);
      Frho.ExtrapolationMethod = 'nearest';
      [X,Y,Z]=meshgrid(linspace(0,1,ceil((grid_rho.ncells)^(1/3))));
      V = Frho(X,Y,Z);
      xslice = 0.5;   
      yslice = [];
      zslice = 0.5;
      slice(X,Y,Z,V,xslice,yslice,zslice,'nearest')
      axis square
      colorbar
		end
	end
