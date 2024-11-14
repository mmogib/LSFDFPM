function [F] = cp11(x)
% Evaluate the  Modified Zhou and Li (2007) function
% x0=+-(0.1,0.1,...) or +-(0.001,0.001,...)'
%x0=((1:1:n)*(1/n+1)-1)'
n = length(x);
%F = zeros(n,1);
%h=1/(n+1);
% last update 6/8/2018
i=2:(n-1);
F(1)=2*x(1)+sin(x(1))-1;
%F(i)=-2*x(i-1)+2*(x(i))+sin(x(i))-1;
F(i)=-x(i-1)+2*(x(i))+sin(x(i))-1;
F(n)=2*x(n)+sin(x(n))-1;
F=F(:);