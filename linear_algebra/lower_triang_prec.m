function y=lower_triang_prec(x,invM,B2,invN)

% y = P x => ( M  0 ) y1 = x1   
%            ( B2 N ) y2   x2 
Np = size(B2,2);
Nr = size(B2,1);
x1 = x(1:Np);
x2 = x(Np+1:Np+Nr);

y=zeros(size(x));
y(1:Np) = invM(x1);
y2 = x2 - B2*y(1:Np);
y(Np+1:Np+Nr) = invN(y2);
