function [F] = pp(x)
% Define the equations as a function
    v1 = 1.0;
    v2 = 1.05;
    v3 = x(1);  % Variable to be solved
    theta1 = 0;
    theta2 =1;
     theta21 = theta2 ;
     
    theta23 = (theta2 -theta3*abs(v1)) ==1; 
%     theta31 = x(1) - x(3);
%     theta32 = x(2) - x(3);
    
    F(1) = 40 * abs(v1) * abs(v2) * sin(theta21) + 20 * abs(v2) * abs(v3) * sin(theta23) - 4;
    F(2) = 20 * abs(v1) * abs(v3) * sin(theta31) + 20 * abs(v3) * abs(v2) * sin(theta32) + 5;
    F(3) = -20 * abs(v1) * abs(v3) * cos(theta31) - 20 * abs(v3) * abs(v2) * cos(theta32) + 40 * abs(v3) * abs(v3) + 4;
F=F(:);

