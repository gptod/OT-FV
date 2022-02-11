function [real_part,imag_part,eigenvectors,ind]= study_eigenvalues(matrix,label,ID)

  if (~exist('ID','var'))
    ID=1;
  end
  
  [eigenvectors,eigenvalues]=eig(full(matrix));
	eigenvalues=diag(eigenvalues);
	abs_eig=abs(eigenvalues);
	[abs_eig,ind] = sort(abs_eig);
	eigenvalues=eigenvalues(ind);
	eigenvectors=eigenvectors(:,ind);
	
  real_part=real(eigenvalues);
  imag_part=imag(eigenvalues);

  fprintf(ID,'real part eigen %s : min=%1.4e max=%1.4e\n', label,min(real_part),max(real_part));
  fprintf(ID,'imag part eigen %s : min=%1.4e max=%1.4e\n', label,min(imag_part),max(imag_part));
  fprintf(ID,'abs       eigen %s : min=%1.4e max=%1.4e\n', label,min(abs_eig),max(abs_eig));
  fprintf(ID,'cond(%s)=%1.4e\n', label,max(abs_eig)/min(abs_eig));
