function [F]=cp9(x)
% Evaluate the Trigexp function  prob 52 Luksan & VIcek 2003
% solution x*=(1,1,...,1)'
% last update 31/10/2018
n = length(x);
%F = zeros(n,1);
i=2:(n-1);
F(1)=3*x(1)^3+2*x(2)-5+sin(x(1)-x(2))*sin(x(1)+x(2));
F(i)=3*x(i).^3+2*x(i+1)-5+sin(x(i)-x(i+1)).*sin(x(i)+x(i+1))+4*x(i)-x(i-1).*exp(x(i-1)-x(i))-3;
F(n)=-x(n-1)*exp(x(n-1)-x(n))+4*x(n)-3;
F=F(:);
%P=max(x,0);