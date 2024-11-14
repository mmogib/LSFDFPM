function [F] = cp17(x)
% Evaluate Problem 4 in Ding  et al. 2017 otimization
% x0=(100,100,...)
n = length(x);
%F = zeros(n,1);
%F = zeros(n,1);
%h=1/(n+1);
% last update 30/10/2018
i=2:(n-1);
F(1)=-4+4*x(1)*((x(1)^2)+x(2)^2); 
F(i)=-4+ 4*x(i).*((x(i-1).^2)+x(i).^2)+4*x(i).*((x(i).^2)+x(i+1).^2);
F(n)=4*x(n)*(x(n-1)^2+x(n)^2);
F=F(:);