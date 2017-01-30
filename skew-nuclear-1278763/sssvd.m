function [U,S,V]=sssvds(A,k,varargin)
% SSSVDS Skew-Symmetric SVD
% Compute the first few singular values and vectors of a skew-symmetric
% matrix.

optsu = struct(varargin{:});
k = k + mod(k,2); % round up to the nearest odd number
p = 2*k;
maxit = 300;

op = false;
if iscell(A)
    Aop = A{1};
    n = A{2};
    op = true;
else
    n = size(A,1);
end    

% get the process started
v0 = randn(n,1);
[V,T] = slanczos(A,v0,p,'full');

for iter=1:maxit
    % compute the Schur fac
end

if nargout<2
end