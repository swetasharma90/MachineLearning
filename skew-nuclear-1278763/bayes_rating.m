function mu = bayes_rating(R,thresh,dim)
% BAYES_RATING Return the bayes rating for items in a set of ratings.
%
% The bayes rating enforces a threshold on the number of required ratings
% in a dataset.  For items that meet that threshold, it scores them as a
% convex combination of the mean rating of the entire dataset 
% and the item rating.  This approach is what IMDB uses for its top 250 list.
%
% Formally, the bayes rating for threshold k is:
% score(i) = mu*k/(ratings(i) + k) + sum_of_ratings(i)/(ratings(i)+k)
% where mu is the grand-mean rating over all user-item pairs
% ratings(i) is the number of ratings on item i
% sum_of_ratings(i) is the total sum of ratings on item(i)
% and k is the rating threshold.
%
% mu = bayes_rating(R,10) returns the bayes rating of all items with at
% least 10 ratings.  
%
% mu = mean_rating(R,10,2) returns the mean non-zero value in the rows of R.
%
% All items with less than 10 ratings (in these examples) will get a score
% of 0.
%
% Example:
%   

% David F. Gleich
% Copyright, Sandia National Laboratories, 2011

% History
% :2011-02-00: Initial coding

grand_mean = mean(nonzeros(R));

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
valid = n>=thresh;

mu(valid) = mu(valid)./(n(valid)+thresh) + ...
    thresh*grand_mean./(n(valid)+thresh);
%mu(~valid) = min(mu(valid))-1;
mu(~valid) = 0;
