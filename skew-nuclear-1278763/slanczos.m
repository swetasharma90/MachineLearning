function [V,T] = slanczos(A,v,k,varargin)
% LANCZOS Run the Skew Lanczos process for k steps.
%
% The skew lanczos process is identical to the Lanczos process but
% customized for a skew-symmetric matrix A.
%
% T = slanczos(A,v,k) returns the (k+1,k) tridiagonal matrix T with the
% skew Lanczos coefficients.  This call does not store the Lanczos vectors
% internally (and therefore, is more efficient than the 'full' call).
%
% [T,R] = slanczos(A,v,k) returns the (k+1,k) tridiagonal matrix T along
% with two vectors in R necessary to restart the process.  This call does 
% not store the Lanczos vectors internally.
%
% [V,T] = slanczos(A,v,k,'full') returns the k+1 set of Lanczos vectors and in
% opposite order.  This call does store the Lanczos vectors internally.
%
% [T,R] = lanczos(A,v,k,T,R) continues a Lanczos factorization for k steps.
% [V,T] = lanczos(A,v,k,T,V) continues a Lanczos factorization for k steps.
% [V,T] = lanczos(A,v,k,T,V,'full') see note 2 below.
%
% Note 1: in the "extend the factorization" calls, the vector v must be
% correctly sized for A, it does not have to be the starting vector for the
% Lanczos process.
%
% Note 2: There is a slight ambiguity in the "update" the factorization
% call.  In particular, if we try to update the factorization after a
% single step [V,T] = lanczos(A,v,1,'full'); [V,T] = lanczos(A,v,1,T,V);
% then the algorithm will erronously detect the second call as [T,R] =
% lanczos(A,v,k,T,R).  This occurs because the file uses size(V,2) == 2 to
% detect the [T,R] input instead of the [T,V] input.  To clarify the
% ambiguity, just specify the 'full' on the end.
%
% The matrix A must be an input that is a valid singleton input to the
% matrixop class or a matrixop class.  See matrixop for more details.
%
% See also LANCZOS

% Copyright David F. Gleich, 2009
% University of British Columbia

% History
% :2009-11-22: Initial coding based on lanczos.m

% A must be an input that creates a matrix operator
Aop = A;

n = length(v);

fullmode = 0;
extend = 0;

% parse the arguments
if isempty(varargin)
    % fall through, just here for simplicity of code
elseif length(varargin) == 1
    % the only option they could have specified is mode...
    fullmode = 1;
elseif length(varargin) == 2
    extend = 1;
    T = varargin{1};
    R = varargin{2};
    if (size(R,2) > 2)
        V = R;
        fullmode = 1;
    end
elseif (length(varargin) == 3)
    if (~strcmpi(varargin{3},'full'))
        error('lanczos:invalidArgument',...
            'invalid third optional argument (it should be ''full'')');
    end
    extend = 1;
    T = varargin{1};
    V = varargin{2};
    fullmode = 1;
else
    error('lanczos:invalidArgument',...
            'Too many input arguments.');
end
    
if ~fullmode
    if (~extend)
        % in this case, we are starting the factorization
        v = v./norm(v);
        b = zeros(k,1);

        b1 = 1;
        v1 = 0;

        ii=0;
        start = 0;
    else
        % in this case, we are extending the factorization
        v = R(:,2);
        v1 = R(:,1);
        
        bstart = -diag(T,-1);

        b = [bstart; zeros(k,1)];
        
        ii = min(size(T));
        b1 = bstart(end);
        start = ii;
    end
    
    while ii < k+start
        ii=ii+1;
        p = Aop*v;
        p = p - b1*v1;
        b(ii) = norm(p);
        b1 = b(ii);
        v1 = v;
        v = -p/b1;
    end;

    bend = [0; b(1:end-1)];
    V = spdiags([-b bend], [-1 1], k+start+1,k+start);
    T = [v1 v];
else
    
    % this mode saves all the vectors
    if (~extend)
        v = v./norm(v);
        b = zeros(k,1);

        V = zeros(n,k+1);

        b1 = 1;
        v1 = 0;

        ii=0;
        start = 0;
    else
        Vstart = V;
        V = [Vstart zeros(n,k)];
        
        % extract all the information from the stored vectors
        v = Vstart(:,end);
        v1 = Vstart(:,end-1);
        
        bstart = -diag(T,-1);
        b = [bstart; zeros(k,1)];
        
        ii = min(size(T));
        b1 = bstart(end);
        start = min(size(T));
    end
    
    while (ii < k + start)
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
end;