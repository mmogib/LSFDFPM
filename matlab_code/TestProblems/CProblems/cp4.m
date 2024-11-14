function [F] = cp4(x)
% Evaluate the Discrete boundary value problem
% for convex constraints 
% last update 24/04/2019
% Problem 7 of RAIRO
n = length(x);
%F = zeros(n,1);
h=1/(n+1);
i=2:(n-1);
%t(i)=i*h;
F(1)=2*x(1)+0.5*h^2*(x(1)+h)^3-x(2);
%F(1)=2*x(1)+0.5*h^2*(x(1)+h)^3;
F(i)=2*x(i)-x(i-1)+x(i+1)+0.5*(h^2)*((x(i)+i'.*h).^3);
F(n)=2*x(n)-x(n-1)+0.5*h^2*(x(n)+(n*h))^3;
F=F(:);
%P=max(x,0);

%Discrete boundry value problem

% @article{cruz2003nonmonotone,
%   title={Nonmonotone spectral methods for large-scale nonlinear systems},
%   author={Cruz, William La and Raydan, Marcos},
%   journal={Optimization Methods and software},
%   volume={18},
%   number={5},
%   pages={583--599},
%   year={2003},
%   publisher={Taylor \& Francis}
% }