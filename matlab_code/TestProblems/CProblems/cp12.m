function [F]=cp12(x)
% Evaluate the discretize bidimensional laplacian operator
% n must be perfect square
% last update 31/10/2018
n = length(x);
r=sqrt(n);
h=1/(r+1);
%F = zeros(n,1);
i=1:(n);
B=(gallery('tridiag',r,-1,4,-1));
A=kron(-eye(r),B);
F(i)=A*x+h^2*(x(i).^3-10);
F=F(:);