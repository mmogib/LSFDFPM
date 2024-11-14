function [F] = cp141(x)
% Pursuit-Evasion problem 21/1/19
n = length(x);
F = zeros(n,1);
i=1:(n);
F(i)=((8)^0.5)*x(i)-1;
F=F(:);
%P=max(x,0); % projection