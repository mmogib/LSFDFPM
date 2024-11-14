% Example 21
function result = F21_DEL(x)
    n = length(x);
    e = ones(n,1);
    A = spdiags([-e 2*e -e],-1:1,n,n);
    f = exp(x) - 1;
    result = A * x + f;
end
