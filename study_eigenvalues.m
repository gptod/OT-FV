function [eigenvalues]= study_eigenvalues(matrix,label,ID)

  if (~exist('ID','var'))
    ID=1;
  end
  
  eigenvalues=eig(matrix);
  real_part=real(eigenvalues);
  imag_part=imag(eigenvalues);
  abs_eig=abs(eigenvalues);
  fprintf(ID,'real part eigen %s : min=%1.4e max=%1.4e\n', label,min(real_part),max(real_part));
  fprintf(ID,'imag part eigen %s : min=%1.4e max=%1.4e\n', label,min(imag_part),max(imag_part));
  fprintf(ID,'abs       eigen %s : min=%1.4e max=%1.4e\n', label,min(abs_eig),max(abs_eig));
  
