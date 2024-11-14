function [F] = cp7(x)
% Evaluate the tridiagonal exponential function 
% x0=(1.5,1.5,...)
n = length(x);
%F = zeros(n,1);
h=1/(n+1);
i=2:(n-1);
F(1)=x(1)-exp(cos(h*(x(1)+x(2))));
F(i)=x(i)-exp(cos(h.*(x(i-1)+x(i)+x(i+1))));
F(n)=x(n)-exp(cos(h*(x(n-1)+x(n))));
F=F(:);



%Tridiagonal function
% Problem 6 in 
% @article{cruz2017spectral,
%   title={A spectral algorithm for large-scale systems of nonlinear monotone equations},
%   author={Cruz, William La},
%   journal={Numerical Algorithms},
%   volume={76},
%   pages={1109--1130},
%   year={2017},
%   publisher={Springer}
% }