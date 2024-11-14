function ls = LSVIIFun(f,bita,ro,sig,varargin)
if nargin < 5, ro = 0.5; end
if nargin < 4, sig = 0.001; end
if nargin < 3, bita = 1; end
ls = @line_search;
    function [alpha,bck] = line_search(x,d,norm_d,varargin)
        % bita = 1; ro = 0.5; sig = 0.001; r=0.1; psi=0.2; gn=1.8;
        % Default parameter values
        if nargin < 3
            aac = matlab.lang.correction.AppendArgumentsCorrection({'f','x','d','norm_d'});
            error(aac, 'MATLAB:notEnoughInputs', 'Not enough input arguments.')
        end
        bck = 0;
        alpha=1;
        xtrial = x + alpha * d;
        Ftrial = f(xtrial);
        Ftd = Ftrial'*d;
        while true
            if bck > 10000
                return;
            end
            if -Ftd >= sig*bita*(norm_d^2)
                return;
            else
                alpha = alpha/2.0;
            end
            xtrial = x + alpha * d;
            Ftrial = f(xtrial);
            Ftd = Ftrial'*d;
            bck = bck + 1;
        end
    end
end