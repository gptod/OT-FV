function y=upper_triang_prec(x,invM,B1T,invN)

% y = P x => ( M B1T ) y1 = x1   
%            ( 0 N   ) y2   x2
Np = size(B1T,1);
Nr = size(B1T,2);
x1 = x(1:Np);
x2 = x(Np+1:Np+Nr);

y2 = invN(x2);
y1 = x1-B1T*y2;
y1 = invM(y1);

y=[y1;y2];
