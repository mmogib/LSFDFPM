function [F] = cp2(x)
% Evaluate the modofied logarithmic function 1 La Cruz et. al 2004
% for solving convex constraints monotone equations
%x0=(1,1,...1)
% to call the function evaluation and projection type [P,F]=feval(f,x)
% to call function evaluation only type [~,F]=feval(f,x)
% last update 26/03/2018
n = length(x);
i=1:(n);
F(i)=log(x(i)+1)-(x(i)/n);
F=F(:);
%P=max(x,0);


%problem is Logarithmic function


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