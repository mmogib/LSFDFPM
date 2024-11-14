function [F] = cp5(x)
% Evaluate the Strictly convex function 1 La Cruz et. al. 2004
% x0=(1/n,2/n,...,1)' (1, 1/2, 1/3,.., 1/n)'
% to call the function evaluation and projection type [P,F]=feval(f,x)
% to call function evaluation only type [~,F]=feval(f,x)
% last update 26/03/2018
n = length(x);
i=1:(n);
F(i)=exp(x(i))-1;
F=F(:);
%P=max(x,0);


% Example 4.1 in
% @article{wang2007projection,
%   title={A projection method for a system of nonlinear monotone equations with convex constraints},
%   author={Wang, Chuanwei and Wang, Yiju and Xu, Chuanliang},
%   journal={Mathematical Methods of Operations Research},
%   volume={66},
%   pages={33--46},
%   year={2007},
%   publisher={Springer}
% }