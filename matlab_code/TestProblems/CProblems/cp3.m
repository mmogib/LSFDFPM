function [F] = cp3(x)
% Evaluate a nonsmooth and monotone function
%
n = length(x);
%F = zeros(n,1);
%h=1/(n+1);
i=1:(n);
F(i)=2*x(i)-sin(x(i));
F=F(:);
%P=max(x,0);



% Problem 2, in  
% @article{zhou2007limited,
%   title={Limited memory BFGS method for nonlinear monotone equations},
%   author={Zhou, Weijun and Li, Donghui},
%   journal={Journal of Computational Mathematics},
%   pages={89--96},
%   year={2007},
%   publisher={JSTOR}
% }