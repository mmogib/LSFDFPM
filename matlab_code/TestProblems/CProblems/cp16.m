function [F] = cp16(x)
% Evaluate Problem Penalty 1
% x0=(100,100,...)
n = length(x);
F = zeros(n,1);
i=1:(n);
c=10^-5;
j=1:(n);
F(i)=2*c*(x(i)-1)+4*(sum(x(j))-1)*x(i);
F=F(:);