function [ITRS, FEVALS,Etime,NORM_F]= runFun(solvernum, num_prob,num_dim, inpointnum)
cons=inpointnum;
zs = num_prob*num_dim*inpointnum;
ITRS = NaN(zs,1);
FEVALS = NaN(zs,1);
Etime = NaN(zs,1);
NORM_F = NaN(zs,1);
tol = 1e-6;maxfev = 2000; maxit = 2000;
counter=1;
% solvers = {'NoLineSearch','LSI','LSII','LSIII','LSV'}
for i=1:num_prob
    for j=1:num_dim
        for l=1:inpointnum
%             k=cons*(j-1)+l+2; 
            probnum=i;
            dimnum=j;
            inpointnum=l;
%             xlrange=xlRange(solvernum,k);
            print_summary(true,false,'solver',solvernum,'func',num2str(i),'dim',num2str(j), 'end_point',num2str(inpointnum));
            switch solvernum
                case 'NoLineSearch'
                   [iterations, func_evals,bck,etime,norm_f] = withoutLS(probnum,dimnum,inpointnum,tol,maxit,maxfev); % Modified FENG method by Auwal and Kumam 2019
%                 case 'LSI'
%                    [iterations, func_evals,bck,etime,norm_f] = LSI(probnum,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case 'LSI'
                   [iterations, func_evals,bck,etime,norm_f] = withLS(probnum,@LSIFun,@SD1Fun,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case 'LSII'
                    [iterations, func_evals,bck,etime,norm_f] = withLS(probnum,@LSIIFun,@SD1Fun,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case 'LSIII'
                    [iterations, func_evals,bck,etime,norm_f] = withLS(probnum,@LSIIIFun,@SD1Fun,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case 'LSIV'
                    [iterations, func_evals,bck,etime,norm_f] = withLS(probnum,@LSIVFun,@SD1Fun,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case 'LSV'
                    [iterations, func_evals,bck,etime,norm_f] = withLS(probnum,@LSVFun,@SD1Fun,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case 'LSVI'
                    [iterations, func_evals,bck,etime,norm_f] = withLS(probnum,@LSVIFun,@SD1Fun,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case 'LSVII'
                    [iterations, func_evals,bck,etime,norm_f] = withLS(probnum,@LSVIIFun,@SD1Fun,dimnum,inpointnum,tol,maxit,maxfev); % PCG 
                case '23'
                    withoutLSN(probnum,dimnum,inpointnum); % PDY Liu and Feng 2018 
                otherwise
                    disp('wrong solver number')
                    break
            end
            if isnan(norm_f)
                iterations=nan; func_evals=nan;bck=nan;etime=nan;
            end
            print_summary(false,true,'itrs',num2str(iterations),'fevals',num2str(func_evals),'bck',num2str(bck), ...
                       'time',num2str(etime),'f_norm',num2str(norm_f));            
            ITRS(counter,1)=iterations;
            FEVALS(counter,1) =func_evals; 
            Etime(counter,1)=etime;
            NORM_F(counter,1)=norm_f;
            counter = counter+1;
        end
    end
end
