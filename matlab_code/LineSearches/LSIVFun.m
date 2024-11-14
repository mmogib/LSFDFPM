function ls = LSIVFun(f,bita,ro,sig,varargin)
if nargin < 5, ro = 0.5; end
if nargin < 4, sig = 0.001; end
if nargin < 3, bita = 1; end
maxitrs=varargin{1};
b = zeros(maxitrs + 1,1);
for i =1:maxitrs+1
    if 1 == i
        b(i) = 1;
    elseif i == 2
        b(i) = 0.5;
    else
        b(i) = 0.5 * sum(b(i-2:i-1));
    end
end
lambda = b(2:end);
ls = @line_search;
    function [alpha,bck] = line_search(x,d,norm_d,varargin)
        % bita = 1; ro = 0.5; sig = 0.001; r=0.1; psi=0.2; gn=1.8;
        % Default parameter values
        k = varargin{1};
       
        if nargin < 3
            aac = matlab.lang.correction.AppendArgumentsCorrection({'f','x','d','norm_d'});
            error(aac, 'MATLAB:notEnoughInputs', 'Not enough input arguments.')
        end
        bck = 0;
        while true
            if bck > 10000
                return;
            end
            alpha = ro^bck;
            Fnew = f(x + bita * alpha * d);
            lhs = -Fnew' * d;
            nFnew =norm(Fnew);
            gg1 = lambda(k) + (1 - lambda(k)) * nFnew;
            rhs = sig * bita * alpha * gg1 * norm_d^2 ;
            if lhs >= rhs
                return;
            end
            bck = bck + 1;
        end
    end
end