function [F] = cp20(x)
% Evaluate the Variational inequality problem Han (2003) A new hybrid generalized 
% proximal point algorithm for variational inequality problems J Global Optim 26
% length of x is n, that is F is from R^n to R^n
% Convex set C=R^n+
% F(x)=D+Mx+q
% D=b*arctan(x),b is a random variable in (0,100) 
% M=A'A+B, A is nxn matrix whose entries
% are randomly generated in (-1,1) and B is a skew-symmetrix matrix
% generated in the same way as A
% q is uniformly distribution in the interval (-500,500)
% rho is constant
% last modified 30/10/2018
n = length(x);
m=500;
%F = zeros(n,1);
%S=rand(n);
A=-1+2.*rand(n);
S=A;
M=(A'*A)+0.5*(S-S');
q=-m+(2*m).*rand(n,1);
b=100*rand(n,1);
i=1:n;
D=b(i).*atan(x(i));
F=D(i,:)+M(i,:)*x+q(i);

F=F(:);
