function vector_portion = get_slice(vector,portion, length_slide, number_of_slide)

				% check dimension
  
	 
  if ( length_slide*number_of_slide ~=  size(vector,1) )
    disp(length_slide)
    disp(number_of_slide)
    
    disp('error in get_slide ')
    return
  end

  
  vector_portion = vector((portion-1)*length_slide+1 : portion*length_slide  );

end
