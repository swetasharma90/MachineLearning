function [A,b] = consmatcomp(varargin)
% CONSMATCOMP Build matrix completion constraints
% [A,b] = consmatcomp(i,j,v,m,n) builds constraints to complete a matrix
% starting from entries given by the triplet arrays i,j,v
% [A,b] = consmatcomp(A0) builds constraints to complete a matrix where A0
% is already a sparse matrix.
% [A,b] = consmatcomp([i,j,v],m,n) builds constraints to complete a matrix
% from entries given by the triplet vectors.  
%
% The third input method should not be used with 3x3 matrices without
% specifying m and n.

% Copyright, David F. Gleich and Lek-Heng Lim
% October 2009

% handle a single argument
if nargin==1
    A0 = varargin{1};
    if size(A0,2) == 3
        % could be a matrix too, check for valid indicies
        try
            A0t = sparse(A0(:,1),A0(:,2),A0(:,3));
        catch 
            % not a valid triplet, must be a matrix
            [i j v] = find(A0);
            A0(:,1) = i; A0(:,2) = j; A0(:,3) = v;
        end
        i = A0(:,1); j = A0(:,2); v = A0(:,3);
        m = max(i); n = max(j);
    else
        % must be a matrix
        [i j v] = find(A0);
        m = size(A0,1); n = size(A0,2);
    end
elseif nargin==2
    % must be [i,j,v],m
    i = A0(:,1); j = A0(:,2); v = A0(:,3);
    m = varargin{2};
    n = max(i);
    if any(i<0) || any(i>m)
        error('invalid dimensions');
    end
elseif nargin==3
    % could be i,j,v or [i,j,v],m,n
    i = varargin{1}; j = varargin{2}; v = varargin{3};
    m = max(i); n = max(j);
    % Figure this out later
elseif nargin==4
    i = varargin{1}; j = varargin{2}; v = varargin{3};
    m = varargin{4}; n = max(j);
elseif nargin==5
    i = varargin{1}; j = varargin{2}; v = varargin{3};
    m = varargin{4}; n = varargin{5};
end

b = v;
A = sparse(1:numel(v), sub2ind([m,n],i,j), 1, numel(v), m*n);





    

    
    
    
                