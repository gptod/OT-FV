function [fileID] = write2td_sequence( fileID, data ,time,section)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (section == 'head')
    fprintf(fileID,'1 %d \n',size(data,1));
end

if (section == 'body')
    fprintf(fileID,'time %9.2e\n',time );
    fprintf(fileID,'%d\n',size(data,1));
    for i=1:size(data,1)
        fprintf(fileID,'%d %16.8e \n',i,data(i));
    end
end
if (section == 'tail')
    fprintf(fileID,'time 1.0e30');
    fclose(fileID);
end 

end


