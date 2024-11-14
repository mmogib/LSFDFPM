function [F] = cp10(x)
% Evaluate the semismooth function Yamashita and Fukushima 1997 Modified
% Newton method for solving semismooth... Math. Program. 76
% length of x is 4, that is F is from R^4 to R^4
% x*=(2,0,1,0)'
% last modified 30/10/2018
%n = length(x);
%F = zeros(n,1);
%A=[1 0 0 0; 0 1 -1 0; 0 1 1 0; 0 0 0 0];
%b=[-10;1;-3;0];
%y=[x(1)^3;x(2)^3;2*x(3)^3;2*x(4)^3];
%F=zeros(n,1);
%i=1:n;
%F=A(i,:)*x+y(i)+b(i);
F(1)=x(1)+x(1)^3-10;
F(2)=x(2)-x(3)+x(2)^3+1;
F(3)=x(2)+x(3)+2*x(3)^3-3;
F(4)=2*x(4)^3;
F=F(:);
