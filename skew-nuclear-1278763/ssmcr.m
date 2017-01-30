function [s,res,mcres] = ssmcr(A,varargin)
% SSMCR Skew-symmetric matrix-completion ranking
%
% s = ssmcr(A) computes a representative set of scores for the items in a
% users-by-items rating matrix A.  This function is a composition of the
% other functions in this package.
%
% [s,res,mcres] = ssmcr(A) returns the:
%   score vector s, the residual score res, and the matrix-completion
%   residual mcres
%
% res, and mcres are not yet returned.  They are coming soon.
%
% ... = ssmcr(A,options_struct) 
% or
% ... = ssmcr(A,'option',value,'option',value...) 
% sets optional arguments.  The struct syntax assumes that 
%
% The options are:
%  skewtype: [ {am} | gm | bc | sb | lo ]
%    The type of skew-symmetric matrix to complete, see ratings2skew.
%  minpair: interger >= 1, default = 1
%    The minimum number of pairwise comparisons required to use an element
%    for matrix completion.
%

options = [];
options.skewtype = 'am';
options.minpair = 1;

if length(varargin) == 1
    optionsu = varargin{1};
elseif length(varargin) > 1    
    optionsu = struct(varargin{:});
else
    optionsu = [];
end
for opt=fieldnames(options)'
    if isfield(optionsu,opt{1}), options.(opt{1}) = optionsu.(opt{1}); end
end

[Y,N] = ratings2skew(A,options.skewtype);
if options.minpair>1, 
    N = spfun(@(x) x>=options.minpair, N);
end
[U,S,V] = svpmc({Y,N},2);
s = U*sum(V,1)';

%mcres = U*S*V