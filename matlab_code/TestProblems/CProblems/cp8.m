function [F] = cp8(x)
% Evaluate the monotone and nonsmooth (at x=1) function
% Yu, Niu & Ma 2013 JIMO
n = length(x);
%F = zeros(n,1);
%h=1/(n+1);
i=1:(n);
F(i)=x(i)-sin(abs(x(i)-1));
F=F(:);
%P={x\in R^n: sum(x)<=n and x>=-1} 


%Nonsmooth function Problem 2

% @article{sun2015modified,
%   title={A modified Hestenes--Stiefel projection method for constrained nonlinear equations and its linear convergence rate},
%   author={Sun, Min and Liu, Jing},
%   journal={Journal of Applied Mathematics and Computing},
%   volume={49},
%   pages={145--156},
%   year={2015},
%   publisher={Springer}
% }