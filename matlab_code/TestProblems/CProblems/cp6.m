function [F] = cp6(x)
% Evaluate the Strictly convex function 2 La Cruz et. al. 2004
% x0=(1/n,2/n,...,1)' (1, 1/2, 1/3,.., 1/n)'
% to call the function evaluation and projection type [P,F]=feval(f,x)
% to call function evaluation only type [~,F]=feval(f,x)
% last update 03/11/2018
n = length(x);
i=1:(n);
F(i)=(i/n)'.*exp(x(i))-1;
F=F(:);
%P=max(x,0);

%Strictly convex function II

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