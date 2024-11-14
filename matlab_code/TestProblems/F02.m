function [F] = F02(x)

n = length(x);
i=1:(n);
F(i)=log(x(i)+1)-(x(i)/n);
F=F(:);


end
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