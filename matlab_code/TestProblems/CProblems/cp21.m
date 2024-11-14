function [F] = cp21(x)
% Evaluate the Variational inequality problem Han (2003) A new hybrid generalized 
% proximal point algorithm for variational inequality problems J Global Optim 26
% length of x is 5, that is F is from R^5 to R^5
% Convex set C={x\in R^5: sum(x)>10, min(x)>=0}
% F(x)=rhoD+Mx+q
% D=arctan(x-2), M is a 5x5 asymmetric +ve definite matrix whose entries
% are randomly generated in (-5,5)
% q is uniformly distribution in the interval (-10,10)
% rho is constant
% last modified 30/10/2018
n = length(x);
ro=1;
m=10;
%F = zeros(n,1);
%S=rand(n);
S=-n+(2*n).*rand(n);
M=0.5*(S-S');
q=-m+(2*m).*rand(n,1);
i=1:n;
D=atan(x(i)-2);
F=ro*D(i,:)+M(i,:)*x+q(i);
F=F(:);
