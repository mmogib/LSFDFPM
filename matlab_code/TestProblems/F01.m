function [F] = F01(x)
n = length(x);
i=1:(n);
F(i)=2*x(i)-sin(x(i));
F=F(:);
end

% Problem 2, in  
% @article{zhou2007limited,
%   title={Limited memory BFGS method for nonlinear monotone equations},
%   author={Zhou, Weijun and Li, Donghui},
%   journal={Journal of Computational Mathematics},
%   pages={89--96},
%   year={2007},
%   publisher={JSTOR}
% }