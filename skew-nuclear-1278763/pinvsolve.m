function X=pinvsolve(A,B,tol)
% PINVSOLVE Solve a system with the pseudo-inverse
% Y = pinvsolve(A,B)
% returns Y = pinv(A)*B


%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.12.4.2 $  $Date: 2004/12/06 16:35:27 $
%  Revised by David F. Gleich

% History
% 2009-12-05: Initial coding based on the mathworks function

if isempty(A)     % quick return
  Y = zeros(size(A',1),size(B,2),class(B));  
  return  
end

[m,n] = size(A);

[U,S,V] = svd(A,0);
if m > 1, s = diag(S);
  elseif m == 1, s = S(1);
  else s = 0;
end
if nargin == 3
  tol = varargin{1};
else
  tol = max(m,n) * eps(max(s));
end
r = sum(s > tol);
s = diag(ones(r,1)./s(1:r));
X = V(:,1:r)*(s*(U(:,1:r)'*B));