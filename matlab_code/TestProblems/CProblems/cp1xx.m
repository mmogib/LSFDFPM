function [F] = cp1(x)
% Evaluate the modified exponential function 2 of La Cruz et. al. 2004
% for convex constraints
% x0=(1/(n^2),1/(n^2),...)' 
% to call the function evaluation and projection type [P,F]=feval(f,x)
% to call function evaluation only type [~,F]=feval(f,x)
% last update 30/03/2018
n = length(x);
F = zeros(n,1);
i=2:(n);
F(1)=exp(x(1))-1;
F(i)=exp(x(i))+x(i)-1;
F=F(:);
%P=max(x,0); % projection


% Strictly convex function I