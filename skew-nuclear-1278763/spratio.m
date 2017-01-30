function [i,j,v]=spratio(A1,A2)
% SPRATIO Compute a triplet form for the ratio of two sparse matrices
% [i,j,v]=spratio(A1,A2) computes the triplet form of A1./A2 
% but only over the non-zeros of A2.
% B = spratio(A1,A2) = A1./A2 over the non-zeros of A2.
%

% Copyright, David F. Gleich 2009
% University of British Columbia

% History
% 2009-11-26: Initial coding
% 2009-12-05: Added sparse return for a single output.
% 2011-02-16: Added better handling of full matrices

m = size(A2,1); n = size(A2,2);

if ~issparse(A2)
    if ~issparse(A1), A1 = full(A1); end
    v = A1(A2~=0)./A2(A2~=0);
    [i,j]=find(A2);
else
    [i,j,v2]=find(A2); % get the full non-zero structure

    % build a map between (i,j) -> non-zero element in v2
    Ai=sparse(i,j,1:length(v2),size(A2,1),size(A2,2));
    a1i = nonzeros((A1~=0).*Ai); % find the indices of non-zeros in A1
    v1 = nonzeros(A1.*spones(A2));
    v = zeros(size(v2));
    v(a1i) = v1;
    v = v./v2;
end

if nargout < 2
    i = sparse(i,j,v,m,n);
end