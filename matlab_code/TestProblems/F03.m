function [F] = F03(x)

n = length(x);
i=1:(n);
F(i)=exp(x(i))-1;
F=F(:);


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
end