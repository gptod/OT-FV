nAi=ncellphi;
denseAi_iplus=zeros(nAi,nAi);
diag_block_invA=sparse_inverse;
if (0) 
      difference=dense_S-C;
				% create blocks
      nAi=ncellphi;
      ctrl_loc=ctrl_solver;
      for i=2:N
			% we need a new seed for agmg, in case is used
	index_agmg=index_agmg+1;
	ctrl_loc=ctrl_solver;
			      % passing index_agmg we will use agmg-i 
	ctrl_loc.init('direct',...
		      controls.ctrl_inner11.tolerance,...
		      controls.ctrl_inner11.itermax,...
		      controls.ctrl_inner11.omega,...
		      controls.ctrl_inner11.verbose,...
		      sprintf('%s A%d',ctrl_inner11.label,i),...
		      index_agmg);
	diag_block_invA.name=sprintf('inverse A%d',i);
      				% define inverse

	fprintf('A %d\n',i)
	iblock=i;
	Ai     =A((iblock-1)*nAi+1 :     iblock*nAi , (iblock-1)*nAi+1 : iblock*nAi);
	iblock=i-1;
	Aiminus=A((iblock-1)*nAi+1 :     iblock*nAi , (iblock-1)*nAi+1 : iblock*nAi);
	
	
	diag_block_invA.init(Ai,ctrl_loc);     
	
	fprintf('%4.2f',0)
	for j=1:nAi
	  
			     %if (mod(i,50)==0)
			     %  fprintf('\b\b\b\b\b\b%4.2f ',i/Nr*100)
			     %  
			     %end
	  
	  iblock=i-1;
%fprintf('%d %d %d \n',(iblock-1)*nAi+1,iblock*nAi, (iblock-1)*nAi+j);
	  column=Aiminus(:,j);
	  
				%size(column)
	  denseAi_iplus(:,j)=diag_block_invA.apply(column);
	end
	fprintf('\n')

	similar=sparse(1:nAi,1:nAi,(spdiags(Aiminus,0)./spdiags(Ai,0))',nAi,nAi);
	eigenvalues=eig(denseAi_iplus-similar);
	real_part=real(eigenvalues);
	imag_part=imag(eigenvalues);
	abs_eig=abs(eigenvalues);
	fprintf('(A_{i+1})^{-1} A_i-I : max/min=%1.4e\n', max(abs_eig)/min(abs_eig))
	fprintf('real part eigen (A_{i+1})^{-1} A_i-I : min=%1.4e max=%1.4e\n', min(real_part),max(real_part))
	fprintf('imag part eigen (A_{i+1})^{-1} A_i-I : min=%1.4e max=%1.4e\n', min(imag_part),max(imag_part))
	fprintf('abs       eigen (A_{i+1})^{-1} A_i-T : min=%1.4e max=%1.4e\n', min(abs_eig),max(abs_eig))

				%figure
				%plot(eigenvalues,'o')

				%figure
				%imagesc(denseAi_iplus)

	if (i==2)
	  jbegin=1;
	  jend  =2*ncellphi;
	elseif(i==N)
	  jbegin=(N-1)*ncellphi+1;
	  jend=Nr;
	else
	  jbegin=(i-3)*ncellphi+1;
	  jend=(i+1)*ncellphi;
	end

	fprintf('jbegin %d jend %d %d \n',jbegin,jend,Nr)
	if (0)
	  fprintf('%4.2f\n',0)
	  for j=jbegin:jend
		    %if (mod(j,100)==0)
		    %  fprintf('\b\b\b\b\b\b%4.2f ',j/Nr*100)	  
		    %end
	    
	    iblock=i-1;
		     %fprintf('%d %d \n',(iblock-1)*nAi+1,iblock*nAi);
	    column=L((iblock-1)*nAi+1:iblock*nAi,j);
	    
				%size(column)
	    difference((iblock-1)*nAi+1:iblock*nAi,j)=diag_block_invA.apply(column);
	  end

	end
	
	
      end

    end
