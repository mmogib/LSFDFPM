function [iterations, func_evals, bck, etime,norm_f] = withoutLS(fnum, dimnum, xnum, tol, maxit, maxfev)
% Derivative-Free Projection Method Without Line Search
% ibrahimkarym@gmail.com

% Default parameter values
% if nargin < 7, maxfev = 2000; end
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


%gamma = 0.2; alpha = 0.5; r=0.1; psi=0.2; gn=1;   %d1
%mu=0.1;

gamma = 0.2; alpha = 0.9; r=0.1; psi=0.2;  %d2
mu=0.1;


% Start timing the algorithm
tic;

% Main loop
while true
    % Check function evaluation limit
%     if any(isnan(x0))
%         fprintf('No convergent %d\n', j);
%         iterations=nan; func_evals=nan;etime =nan; norm_f=nan;
%         return;
%     end
    if norm_F0 <= tol 
        %fprintf('Solution (with x) found at iteration %d\n', j);
        iterations=j; func_evals= FEV + bck;etime =toc; norm_f=norm_F0;
        return;
    end


    if j >= maxit
        %disp('Maximum number of iterations reached');
        iterations=nan; func_evals=nan;etime =nan; norm_f=nan;
        return;
    end
    
%     if FEV >= maxfev
%         disp('Maximum number of function evaluations reached');
%         iterations=j; func_evals=FEV + bck;etime =toc; norm_f=nan;
%         return;
%     end


    % Compute new point and evaluate function
    z = x0 + alpha * d0;
    Fz = f(z);
    FEV = FEV + 1;
    norm_Fz=norm(Fz);
    if  norm_Fz<= tol
        %fprintf('Solution (with z) found at iteration %d\n', j);
        iterations=j; func_evals=FEV + bck;etime =toc; norm_f=norm_Fz;
        return;
    else
        x1 = project_Hj(x0, mu, F0, Fz, z);
    end

    % Step 4: Update the step size parameter mu_j+1
    if norm(F0 - Fz) > 0
        mu = min(gamma * norm(x0 - z) / norm(F0 - Fz), mu);
    end



    % Update variables for next iteration
    F1 = f(x1);
    FEV = FEV + 1;

    y = F1-F0;
    s = x1 -x0 + r*y;
    sk1yk1 = s'*y;
    yk1yk1 =y'*y;
    normY =sqrt(yk1yk1);
    Fxkyk1 = F1'*y;
    Fxkdk1 = F1'*d0;
    vk = max(psi*nrmd*normY,norm_F0^2);
    alpI = sk1yk1/yk1yk1;
    alpII = Fxkdk1/vk;
    betak =Fxkyk1/vk;
    d0 = - alpI*F1 + betak*d0 - alpII*y;

    % Update state
    x0 = x1;
    F0 = F1;
    norm_F0 = norm(F0);
    nrmd = norm(d0);
    j = j + 1;
end


end




function proj = project_Hj(x0, mu, F0, Fz, z)
% Computes the projection of y = x - mu * Psi(zj) onto the hyperplane Hj

% Calculate the normal vector to the hyperplane
v = x0 - mu * F0 - z;

% Calculate the point to be projected
y = x0 - mu * Fz;

% Compute the projection of y onto the hyperplane
proj = y - ((y - z)' * v) / (v' * v) * v;
end