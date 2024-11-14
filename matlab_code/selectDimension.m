function dim = selectDimension(dimnum)
    switch dimnum
        case 1, dim = 500;
        case 2, dim = 15000;
        case 3, dim = 75000;
        case 4, dim = 150000;
        otherwise, dim = dimnum;
    end
end

