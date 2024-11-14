function [d] = SD1Fun(d,norm_d,F0,norm_F0,F1,x0,x1,r,psi)
y = F1-F0;
s = x1 -x0 + r*y;
sk1yk1 = s'*y;
yk1yk1 =y'*y;
normY =sqrt(yk1yk1);
Fxkyk1 = F1'*y;
Fxkdk1 = F1'*d;
vk = max(psi*norm_d*normY,norm_F0^2);
alpI = sk1yk1/yk1yk1;
alpII = Fxkdk1/vk;
betak =Fxkyk1/vk;
d = - alpI*F1 + betak*d - alpII*y;
end

