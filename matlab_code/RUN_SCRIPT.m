addpath(genpath("../matlab_code"));
clear global;
problems = 14;
dims = 4;
points = 9;
% PP = [3;7;11;15;19];% 11:;LSII,15::LSIII,19::LSV
% PP = [3;7];% 11:;LSII,15::LSIII,19::LSV
% PPname={'NoLineSearch','LSI','LSII','LSIII','LSV'};
PPname={'NoLineSearch','LSI','LSII', 'LSIII','LSIV', 'LSV', 'LSVI', 'LSVII'};

solvers = numel(PPname);
columns = 4;
zs = problems*dims*points;

M_Itrs = NaN(zs,solvers);
M_Fev = NaN(zs,solvers);
M_Time = NaN(zs,solvers);
M_Fnorm = NaN(zs,solvers);

T = NaN(zs,columns*solvers); 
problemColm = reshape(repelem(arrayfun(@(i) testProblem(i),1:problems),dims*points),[],1);
pointsColm = repmat(reshape(repelem(arrayfun(@(i) sprintf("x%s",num2str(i)),1:points),dims),dims*points,1),problems,1);
dimColm = repmat(repmat(reshape(arrayfun(@(i) selectDimension(i),1:dims),dims,1),points,1),problems,1);
for i = 1:solvers
    [iterations, func_evals,etime,norm_f] = runFun(PPname{i},problems,dims,points);
    M_Itrs(:,i)=iterations(:);
    T(:,4*(i-1)+1)=iterations(:);
    M_Fev(:,i)=func_evals(:);
    T(:,4*(i-1)+2)=func_evals(:);
    M_Time(:,i)=etime(:);
    T(:,4*(i-1)+3)=etime(:);
    M_Fnorm(:,i)=norm_f(:);
    T(:,4*(i-1)+4)=norm_f(:);
end
currentFilePath = mfilename('fullpath');
[currentFolder, ~, ~] = fileparts(currentFilePath);
Filename = sprintf('%s/results/DFPMWLS_NEW%s.xlsx',currentFolder, datestr(now,'mm-dd-yyyy-HH-MM'));
names = cell(1,columns*solvers);

name_counter=1;
col_names = {'Iterations','FunctionEvals','Time','NormF'};
for i=1:solvers
    for j =1:columns
        names{name_counter}=sprintf("solver%s%s",PPname{i},col_names{j});
        name_counter=name_counter+1;
    end
end
names = cellstr(names);
T1 = table(problemColm,pointsColm,dimColm,'VariableNames',["Problem","InitialPoint","Dimension"]);
T2 =array2table(T,"VariableNames",names);
Ttable = [T1 T2];
writetable(Ttable,Filename);

title_names = {'Iterations','Function Evaluations','Time'};
M = {M_Itrs,M_Fev,M_Time};
for k=1:numel(title_names)
    newperf(M{k}, 5, PPname{:},title_names{k});
end