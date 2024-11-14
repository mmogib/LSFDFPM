function [iterations, func_evals,bck,etime,norm_f] = withLS(fnum,linesearch,sd,dimnum,xnum,tol,maxit,maxfev)
% Derivative-Free Projection Method
% ibrahimkarym@gmail.com

% Default parameter values
if nargin < 7, maxfev = 2000; end
if nargin < 6, maxit = 2000; end
if nargin < 5, tol = 1e-6; end
% if nargin < 4, xlrange = []; end
if nargin < 3, xnum = 1; end
if nargin < 2, dimnum = 1; end
if nargin < 1, fnum = 1; end

% Determine problem dimension
dim = selectDimension(dimnum);

% Select the function based on `fnum`
f = testProblem(fnum);
f = str2func(f);


% Select the initial guess based on `xnum`
x0 = getInitialPoint(xnum, dim);

% Initialization
j = 0; FEV = 0; bck = 0;


F0 = f(x0);  % Evaluate the function at initial guess
FEV = FEV + 1;

% Compute initial values
norm_F0 = norm(F0);
d0 = -F0;
nrmd = norm(d0);

% Line search parameters
bita = 1; ro = 0.5; sig = 0.001; r=0.1; psi=0.2; gn=1.8;
% linesearch function 
lsfun=linesearch(f,bita,ro,sig,maxit);

% Start timing the algorithm
tic;

% Main loop
while true

    if norm_F0 <= tol
        %fprintf('Solution (with x) found at iteration %d\n', j);
        iterations=j; func_evals= FEV;etime =toc; norm_f=norm_F0;
        return;
    end


    if j >= maxit
        %disp('Maximum number of iterations reached');
        iterations=nan; func_evals=nan;etime =nan; norm_f=nan;
        return;
    end

    if FEV >= maxfev
        %disp('Maximum number of function evaluations reached');
        iterations=nan; func_evals=nan;etime =nan; norm_f=nan;
        return;
    end

   

    % Line search
    
    [rom,bck] = lsfun(x0,d0,nrmd,j+1);
    FEV =FEV + bck;


    % Compute new point and evaluate function
    z = x0 + rom * d0;
    Fz = f(z);
    FEV = FEV + 1;
    norm_Fz=norm(Fz);
    if  norm_Fz<= tol
        iterations=j; func_evals=FEV;etime =toc; norm_f=norm_Fz;
        return;
    else
        % Update extrapolated point
        zetak = (Fz' * (x0 - z)) / norm(Fz)^2;
        x1 = x0 - gn*zetak * Fz;

    end

    % Update variables for next iteration
    F1 = f(x1);
    FEV = FEV + 1;

    % Compute new direction
    d0 = sd(d0,nrmd,F0,norm_F0,F1,x0,x1,r,psi);

    % Update state
    j = j + 1;
    x0 = x1;
    F0 = F1;
    norm_F0 = norm(F0);
    nrmd = norm(d0);
end

end