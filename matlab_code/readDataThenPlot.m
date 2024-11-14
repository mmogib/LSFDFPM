addpath("../matlab")
Part1 = readtable(".\Results\DFPMWLS_NEW11-12-2024-21-59.xlsx");
% Part2 = readtable("code-withoutlinesearch\withoutLS\results\DFPMWLS_NEW11-10-2024-08-33.xlsx");
PPname={'NoLineSearch','LSI','LSII', 'LSIII','LSIV', 'LSV', 'LSVI', 'LSVII'};
title_names = {'Iterations','Function Evaluations','Time'};
NoLSItrs = Part1.solverNoLineSearchIterations;
LS1Itrs = Part1.solverLSIIterations;
LS2Itrs = Part1.solverLSIIIterations;
LS3Itrs = Part1.solverLSIIIIterations;
LS4Itrs = Part1.solverLSIVIterations;
LS5Itrs = Part1.solverLSVIterations;
LS6Itrs = Part1.solverLSVIIterations;
LS7Itrs = Part1.solverLSVIIIterations;


NoLSFevals = Part1.solverNoLineSearchFunctionEvals;
LS1Fevals = Part1.solverLSIFunctionEvals;
LS2Fevals = Part1.solverLSIIFunctionEvals;
LS3Fevals = Part1.solverLSIIIFunctionEvals;
LS4Fevals = Part1.solverLSIVFunctionEvals;
LS5Fevals = Part1.solverLSVFunctionEvals;
LS6Fevals = Part1.solverLSVIFunctionEvals;
LS7Fevals = Part1.solverLSVIIFunctionEvals;

vNoLSTime = Part1.solverNoLineSearchTime;
LS1Time = Part1.solverLSITime;
LS2Time = Part1.solverLSIITime;
LS3Time = Part1.solverLSIIITime;
LS4Time = Part1.solverLSIVTime;
LS5Time = Part1.solverLSVTime;
LS6Time = Part1.solverLSVITime;
LS7Time = Part1.solverLSVIITime;


M_Itrs = cat(2,NoLSItrs,LS1Itrs,LS2Itrs,LS3Itrs,LS4Itrs,LS5Itrs,LS6Itrs,LS7Itrs);
M_Fev  = cat(2,NoLSFevals,LS1Fevals,LS2Fevals,LS3Fevals,LS4Fevals,LS5Fevals,LS6Fevals,LS7Fevals);
M_Time = cat(2,NoLSTime,LS1Time,LS2Time,LS3Time,LS4Time,LS5Time,LS6Time,LS7Time);
M = {M_Itrs,M_Fev,M_Time};
for k=1:numel(title_names)
    newperf(M{k}, 5, PPname{:},title_names{k});
end
