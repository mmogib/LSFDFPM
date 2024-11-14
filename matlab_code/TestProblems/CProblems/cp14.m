function [F] = cp14(x)
% HEQ  Chandrasekhar H-equation residual
%      Jacobian uses precomputed data for fast evaluation
%
%      [F] = cp8(X) returns the nonlinear residual
%      F 
%
% Be sure and store the correct data in the global array A_heq.
% % last update 7/7/2018, 30/10/2018
global A_heq;
n=length(x);
%c=0.9;
c=0.99;
%c=0.9999;
mu=1:n; mu=(mu-.5)/n; mu=mu';
cc=.5*c/n;
A_heq=ones(n,1)*mu'; A_heq=cc*A_heq'./(A_heq+A_heq');
h=ones(n,1)-(A_heq*x);
ph=ones(n,1)./h;
F=x-ph;

