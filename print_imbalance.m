function [out] = print_imbalance(v,dimblock)
  for i=1:size(v,1)/dimblock
    fprintf('%0.2d - imb=%+1.1e norm=%1.1e\n',...
						i,sum(v(1+(i-1)*dimblock:i*dimblock)),norm(v(1+(i-1)*dimblock:i*dimblock)))
  end
end  
