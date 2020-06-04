function [L,U,D] = effChol(kernel,Xdata,xnew,L,U,D)
%
% Example 1. 
% dist = @(x,y) abs(x - y');
% kernel =@(dis) (1 + dis).*exp(-dis);
% Xdata = [0:0.1:0.6 0.8:0.1:1]';
% xnew = 0.7;
% K = kernel(dist(Xdata,Xdata));
% R = chol(K);
% D = diag(diag(R).^2);
% L = R'*inv(diag(diag(R)));
% U = inv(L');
% [L,U,D] = effChol(kernel,Xdata,xnew,L,U,D);
dist = @(x,y) abs(x - y');
n = length(Xdata);
ltemp = U'*kernel(dist(Xdata,xnew));
Ltemp = D\ltemp;
dtemp = kernel(dist(xnew,xnew))-ltemp'*Ltemp;
Utemp = -U*Ltemp;
L(n+1,1:n) = Ltemp;
L(n+1,n+1) = 1;
D(n+1,n+1) = dtemp;
U (:,n+1) = Utemp;
U(n+1,n+1) = 1;
