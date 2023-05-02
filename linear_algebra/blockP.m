function y=blockP(x,tol,A,B1,B2,C)

% functin that applies the preconditioner P,
% i.e. solves y=P\x, where P is the right block preconditioner
% P = [S B1;
%      0 C];
% with S=A-B1*C^{-1}*B2

S = A-B1*(C\B2);
Np = size(A,1);
Nr = size(C,1);
x1 = x(1:Np);
x2 = x(Np+1:Np+Nr);

y2 = C\x2;
y1 = B1*y2-x1;

%y1 = -S\y1;
[y1,flag,~,~] = agmg(-S,y1,0,1e-3,3,-1,[],0);
if flag~=0 && flag~=1
    disp('error in the precodintioning')
    fprintf('%7s %1i \n','flag1 = ',flag)
    return
end

y=[y1;y2];



% % left block preconditioner
% % P = [S 0;
% %      B2 C];
% % where S=A-B1*C^{-1}*B2

% S = A-B1*(C\B2);
% Np = size(A,1);
% Nr = size(C,1);
% x1 = x(1:Np);
% x2 = x(Np+1:Np+Nr);
% 
% %y1 = -S\x1;
% [y1,flag,~,~] = agmg(-S,x1,0,1e-3,3,-1,[],0);
% if flag~=0 && flag~=1
%     disp('error in the precodintioning')
%     fprintf('%7s %1i \n','flag1 = ',flag)
%     return
% end
% 
% y2 = C\(B2*y1+x2);
% y1 = -y1;
% 
% y=[y1;y2];