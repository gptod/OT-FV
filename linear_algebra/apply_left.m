function [y] = apply_left(x,invA,B2,Np,Nr)
	y=x;
	y(Np+1:Np+Nr)=y(Np+1:Np+Nr)-B2(invA(y(1:Np)));
end
