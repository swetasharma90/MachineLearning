function [AtA,r]=ssr2(Y,completeflag)
% SSR2 Minimize a Frobenius fit over rank-2 skew-symmetric matrices
%
% [AtA,b]=ssr2(Y) returns the normal equations for the least squares
% problem
%
%    minimize ||Y - (ge^T - eg^T)||_F
%       g
%
% which corresponds with the minimum over rank 2 skew-symmetric matrices.
% By default, only the non-zero values of Y are used.  However,
% [AtA,b]=ssr2(Y,1) changes this behavior to use every single comparison
% with the items of Y.  In the latter case, any 0 values in Y are used
% within the fit.  
%
% The matrix AtA is ALWAYS singular.  Thus, we still need to treat this
% with an appropriate solver.  Alternatively, deleting any non-zero row and
% column will make the matrix non-singular.

% David F. Gleich, Copyright 2009
% University of British Columbia

% History
% :2009-11-22: Initial coding
% :2009-11-26: Switched to solving with AtA instead of A
% :2009-12-06: Made interpretation of Y(i,j) consistent: 
%              Y(i,j) = g(i) - g(j)

if ~exist('completeflag','var'), completeflag=0; end

if completeflag   
    n = size(Y,1);
    AtA = ones(n,n);
    AtA = diag(n*ones(n,1)) - AtA;
    r = zeros(n,1);
    % TODO vectorize this loop
    for i=1:n
        for j=i+1:n
            if i==j, continue; end
            r(i) = r(i) + Y(i,j);
            r(j) = r(j) - Y(i,j);
        end
    end
else
    if issparse(Y)
        n = size(Y,1);
        [i,j,v] = find(Y);
    elseif iscell(Y) && length(Y)==5
        [i,j,v,n,n] = deal(Y{:});
    end
        
    AtA = sparse(i,j,-1,n,n);
    AtA = diag(-sum(AtA)) + AtA;
    r = accumarray(i,v);
    r = r - accumarray(j,v);
    r = r/2; 
end