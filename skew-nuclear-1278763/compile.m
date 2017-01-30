function compile
mex tripletmult.c
mex indexprojfact.c
if isunix
    mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99" ratings2skew_mex.c
else
    mex -largeArrayDims ratings2skew_mex.c
end