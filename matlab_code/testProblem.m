function f = testProblem(fnum)
    problems = ["F01", "F02", "F03", "F04", "F05", "F06", ...
                "F07", "F08", "F09", "F10", "F11", "F12", ...
                "F13", "F14"];
    if fnum <= numel(problems)
        f = problems(fnum);
    else
        f = 'fnum';
    end
end
