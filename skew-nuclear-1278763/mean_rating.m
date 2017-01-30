function mu = mean_ratings(R,dim)
% MEAN_RATING Return the mean rating amongst non-zero values only.
%
% mu = mean_rating(R) returns the mean non-zero value in the columns of R.
%
% mu = mean_rating(R,2) returns the mean non-zero value in the rows of R.
%
% Example:
%   

% David F. Gleich
% Copyright, Sandia National Laboratories, 2011

% History
% :2011-02-04: Initial coding

if ~exist('dim','var') || isempty(dim), dim=1; end
[i,j,v] = find(R);
if dim==1
    ind = j;
    dim = 2;
else
    ind = i;
    dim = 1;
end

mu = accumarray(ind,v,[size(R,dim),1]);
n = accumarray(ind,1,[size(R,dim),1]);
mu = mu./n;
