function out=prec_times_matrix(prec,matrix)

  N=size(matrix,1);
  M=size(matrix,2);
  out=zeros(N,M);

  
  %fprintf('\b\b\b\b\b\b%1.2e ',0.0)
  for i=1:M
    if (mod(i,50)==0)
      %fprintf('\b\b\b\b\b\b\b\b\b%1.2e ',i/M*100)
    end    
    out(:,i)=prec(full(matrix(:,i)));
  end
  %fprintf('\n')
end
