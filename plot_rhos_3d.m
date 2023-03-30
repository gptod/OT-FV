function [] = plot_rhos_3d(grid_rho,rho_in,rho,rho_f,dsp)

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
