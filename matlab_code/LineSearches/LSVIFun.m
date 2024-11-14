function ls = LSVIFun(f,bita,ro,sig,varargin)
if nargin < 5, ro = 0.5; end
if nargin < 4, sig = 0.001; end
if nargin < 3, bita = 1; end
eta = 0.001;
xi = 0.6;
ls = @line_search;
    function [alpha,bck] = line_search(x,d,norm_d,varargin)
        % bita = 1; ro = 0.5; sig = 0.001; r=0.1; psi=0.2; gn=1.8;
        % Default parameter values
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
            P_eta_xi = min(xi, max(nFnew, eta));
            rhs = sig * bita * alpha * P_eta_xi * norm_d^2 ;
            if lhs >= rhs
                return;
            end
            bck = bck + 1;
        end
    end
end