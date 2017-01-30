function [r,nv,i,j,v] = score_residual(A,s)
% SCORE_RESIDUAL Output the residual for a particular score vector
%
% Given a partially filled skew-symmetric matrix Y, the score residual
% is ||Y - (se^T - es^T)||_F (i.e. the frobenius norm), but _only_ over the
% filled/known entries of Y.
% 
% There are four ways to provide the matrix to this function:
% r = score_residual(Y,s); 
% r = score_residual({i,j,v},s);
% r = score_residual({i,j,v,m,n},s); (note that m==n for this call)
% r = score_residual({N,D},s); use non-zero structure of D with valus N./D
%   if D is a logical matrix, this allows zero values in N to be respected
%
% The output is the Frobenius norm of the difference.
%
% [r,nY,i,j,rvec] = score_residual(...) also returns the overall score 
% difference norm (nY = norm(Y,'fro')) along with individual differences and 
% the non-zero entries where they occurred.
%
% This function reports the entire residual, including duplicates if they
% are specified in the matrix.
%
% Example:
%   A = 

% Copyright David F. Gleich, 2011

% History
% :2011-02-16: Initial coding


if issparse(A)
    [i,j,v]=find(A);
    [m,n]=size(A);
    i = uint(i-1);
    j = uint(j-1);
elseif iscell(A)
    if length(A)==2 % compute from A1./A2 but get non-zero structure right
        [m,n]=size(A{2});
        [i,j,v]=spratio(A{1},A{2});
        i=uint32(i-1);
        j=uint32(j-1);
    elseif length(A)==3
        i=A{1}; j=A{2}; v=A{3};
        m = max(i);
        n = max(j);
    elseif length(A)==5
        i=A{1}; j=A{2}; v=A{3};
        m=A{4}; n=A{5};
    end
end
if ~exist('m','var'), error('score_residual:invalidArgument',...
        'please provide a valid matrix input');
end
if m~=n, error('score_residual:invalidArgument', ...
        'the input matrix must be square');
end

nv = norm(v);
for nzi=1:length(i)
    v(nzi) = v(nzi)-(s(i(nzi)+1) - s(j(nzi)+1));
end
r = norm(v)/nv;