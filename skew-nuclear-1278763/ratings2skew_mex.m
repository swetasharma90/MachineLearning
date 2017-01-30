function vargout=ratings2skew_mex(varargin)
% RATINGS2SKEW_MEX

if isunix
    mex -largeArrayDims CFLAGS="\$CFLAGS -std=c99" ratings2skew_mex.c
end

vargout{:} = ratings2skew_mex(varargin{:});