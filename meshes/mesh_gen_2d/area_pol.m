function [area] = area_pol(Np,coord)

%compute the area of the polygon
vertices = [coord; coord(1,:)];
sum1=0; sum2=0;
for i=1:Np
    sum1 = sum1 + vertices(i,1)*vertices(i+1,2);
    sum2 = sum2 + vertices(i,2)*vertices(i+1,1);
end
area = (sum1-sum2)/2;
area = abs(area);

end

