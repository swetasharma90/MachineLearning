function [U S V t err] = svpss (Omega, b, m, n, k, params)
%Return the result of solving min \| b - X(Omega) \|_F^2 , s.t rank(X) <= k
% AND X is skew-symmetric.
%
%Inputs:
%   Omega - indices of known entries.
%   b - values of known entries
%   k - approximate rank. We only work with matrices of rank up to rank k.
%   params - parameter structure containing the following fields:
%       tol,vtol,mxitr - stopping criteria
%           defaults = 1e-2 (tol), 1e-3 (vtol), 100 (mxitr)
%       verbosity - verbosity level
%           default = 1
%       eta - step size
%           default = 3/(4p) where p=fraction of entries available
%Outputs:
%   U - Left singular vectors of the optimal matrix X
%   S - Singular values of X
%   V - Right singular vectors of X
%   t - Number of iterations required
%   err - RMSE over the known values
%
%Reference:
% David F. Gleich and Lek-Heng Lim...
% ...
%
% Raghu Meka, Prateek Jain, Inderjit S. Dhillon. 
% "Guaranteed Rank Minimization via Singular Value Projection"
% Arxiv:0909.5457
%Written by: Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka (raghu@cs.utexas.edu)
%Last updated: October, 2009
% Modified by: David Gleich
% November 2009

error(nargchk(5,6,nargin,'struct'));
if ~isreal(m) || ~isreal(n) || ~isreal(k)
    error('svp:invalidInputs','m, n, and k should all be real integers');
end
if length(Omega) ~= length(b)
    error('svp:invalidDimension','Omega and b should have equal row lengths');
end
if any(Omega)<0
    error('svp:invalidInputs','Omega must be a positive integer matrix');
end
if size(Omega,2)==2
    % Omega is pairs of indices
    if any(Omega(:,1)>m) || any(Omega(:,2))>n
        error('svp:invalidInputs','Omega has invalid indicies (> m or > n)');
    end
else
    if any(Omega>m*n)
        error('svp:invalidInputs','Omega has invalid indicies (> m*n)');
    end
end

p = length(Omega)/(m*n); %sampling density

if ~exist('params','var') || isempty(params), params=struct(); end
params = SetDefaultParams(params,p);
tol=params.tol;
vtol=params.vtol;
mxitr=params.mxitr;
verbose = params.verbosity;
eta=params.eta; %Step Size

% Initialize X=U*S*V'=0
U = zeros(m,k);
S = zeros(k,1);
V = zeros(n,k);

% Compute indices for the known entries
if size(Omega,2) == 2
    I = Omega(:,1);
    J = Omega(:,2);
    Omega = sub2ind([m,n],I,J);
else     
    [I, J] = ind2sub([m,n], Omega);
end

oerr=+Inf;
t = 1;

while  t <= mxitr

    % Compute value of current iterate X=U*S*V' over the known entries
    X_Omega=compute_Xss_Omega(U*diag(S), V, Omega);

    % Compute RMSE and break if RMSE is small enough or not enough progress
    err=norm(X_Omega - b,'fro')/sqrt(length(b));

    if verbose
        fprintf('  %4i  %8.2e  %8.2e\n', t, err, abs(oerr-err));
    end
    if err < tol || abs(oerr - err) < vtol
        break
    end

    oerr = err;

    % Compute a step in direction of gradient, i.e., Y=eta P_omega(M-X)
    y=eta*(b - X_Omega);

    % Take a step in the direction of gradient
    % and compute projection over top k singular vectors
    % i.e, U*S*V'=P_k(U*S*V'+Y)
    [U S V] = fastsvd(U, S, V, y, I, J, k);
    
    t = t+1;

end
fprintf('\n');

function [U,S,V] = fastsvd(U,S,V,y,I,J,k)
%Find the singular values of U S V' + Y
%Y is a sparse matrix
m = size(U,1);
n = size(V,1);
Afunc = @(x) 1/2*((U*(S.*(V'*x)) - V*(S.*(U'*x))) + compOmegaYx(m,n,y,I,J,x));
Atfunc = @(x) 1/2*((U*(S.*(V'*x)) - V*(S.*(U'*x))) + compOmegaYtx(m,n,y,I,J,x));

[U,Sigma,V] = lansvd(Afunc,Atfunc, m,n,k,'L');
S = diag(Sigma);
