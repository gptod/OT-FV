function y=triang_prec(x,invA,B1T,invC)

Np = size(B1T,1);
Nr = size(B1T,2);
x1 = x(1:Np);
x2 = x(Np+1:Np+Nr);

y2 = invC(x(Np+1:Np+Nr)),;
y1 = B1T*y2-x1;
y1 = invA(y1);

y=[y1;y2];
