function [out] = print_imbalance(v,dimblock)
  for i=1:size(v,1)/dimblock
    fprintf('%d - imb=%1.4e norm=%1.4e\n',i,sum(v(1+(i-1)*dimblock:i*dimblock)),norm(v(1+(i-1)*dimblock:i*dimblock)))
  end
end  
