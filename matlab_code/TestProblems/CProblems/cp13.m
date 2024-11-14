function [F]=cp13(x,Extra)
% Evaluate the discretize of nondiffrentiable Dirichlet problem
% n must be perfect square
n = length(x);
r=sqrt(n);
h=1/(r+1);
%F = zeros(n,1);
i=1:(n);
B=(gallery('tridiag',r,-1,4,-1));
A=kron(-eye(r),B);
F(i)=A*x-h^2*((max(x(i)-1,0.5*x(i)-0.5))+ones(n,1));
F=F(:);