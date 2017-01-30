function [T,V] = slanczos_simple(A,v,k)
% SLANCZOS_SIMPLE Run the Skew Lanczos process for k steps.
%
% The skew lanczos process is identical to the Lanczos process but
% customized for a skew-symmetric matrix A.
%
% [T,V] = slanczos_simple(A,v,k) returns the (k+1,k) tridiagonal matrix T along
% with two vectors in R necessary to restart the process.  This call does 
% not store the Lanczos vectors internally.
%
% See also LANCZOS

% Copyright David F. Gleich, 2009
% University of British Columbia

% History
% :2009-11-22: Initial coding based on lanczos.m

% A must be an input that creates a matrix operator
Aop = A;

n = length(v);


% this mode saves all the vectors
v = v./norm(v);
b = zeros(k,1);

V = zeros(n,k+1);

b1 = 1;
v1 = 0;

ii=0;
start = 0;
    
while ii < k + start
    ii=ii+1;
    V(:,ii) = v;
    p = Aop*v;
    p = p - b1*v1;
    b(ii) = norm(p);
    b1 = b(ii);
    v1 = v;
    v = -p/b1;
end
V(:,end) = v;
    
bend = [0; b(1:end-1)];
T = spdiags([-b bend], [-1 1], k+1+start,k+start);





