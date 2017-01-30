function [R,theta,beta]=generate_ratings(nusers,nitems,ratprob,err)
% GENERATE_RATINGS Generate synthetic ratings
%
% R = generate_ratings(nusers,nitems,missing,err) produces a nusers-by-nitems
% rating matrix based on a item response theory model:
%   R(i,j) = rating(alpha(i) + beta(i)*theta(j) + err*eps(i,j))
% where alpha(i) is the center rating for user i, beta(i) is the rating
% sensitivity of user i, and theta(j) is the intrinsic score of item j.
% The eps(i,j) is an error models.  For this function:
%   alpha(i) ~ N(3,1)
%   beta(i) ~ N(0.5,0.5)
%   theta(i) = N(0.1,1)
%   eps(i,j) = N(0,1)
% and rating(score) is a levels function with the following cutpoints:
%  R = 1 if       score <= 1.5
%  R = 2 if 1.5 < score <= 2.5
%  R = 3 if 2.5 < score <= 3.5
%  R = 4 if 3.5 < score <= 4.5
%  R = 5 if 4.5 < score
%
% Also, missing controls the number of missing ratings.  The missing rating
% model biases missing entries to the ratings 1 and 2.  Rating 2 is more
% likely to appear than rating 1, and ratings 3-5 are equally likely to
% appear.  If missing > 1, then it is interpreted as an approximate number
% of movies to rate.  If missing < 1, then it is interpreted as a
% probability of a rating appearing.  
%
% [R,theta,beta,alpha] = generate_ratings(...) also returns the current set
% of scores.
% 
% Example:
%   [R,t] = generate_ratings(10,4,2); % generate ratings 
%    

% Copyright, Sandia Corporation, 2011
% David F. Gleich

% History
% :2011-02-04: Initial coding


if ratprob>1
    ratprob = ratprob/nitems;
end

if ~exist('err','var') || isempty(err), err = 1; end

% draw a center for each user's rating
alpha = randn(nusers,1)+3;
% draw a sensitivity for each user
beta = sqrt(0.5)*randn(nusers,1)+0.5;
%beta = zeros(nusers,1);
% draw a score for each movie
theta = randn(1,nitems)+0.1;

Y = alpha*ones(1,nitems) + beta*theta + err*randn(nusers,nitems);
minY = min(Y(:));
maxY = max(Y(:));
R = ordinal(Y,num2str([1 2 3 4 5]'),[],[minY-1,1.5,2.5,3.5,4.5,maxY]);
R = double(R);
R = R.*((1./(1+2*(R>=3)+(R>=2)).*rand(nusers,nitems))<ratprob);
% = R
%R = Y;