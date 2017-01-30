function top_items(items,scores,topk)
% TOP_ITEMS Report on the best items given a list of scores.
%
% top_items(itemnames, scores) reports the top 15 item names based on the
% largest scores.
%
% top_items(itemnames, scores, k) reports the top k item names.
%
% Example:
%   load movielens_comedy
%   mu = mean_rating(A);
%   top_items(movies, mu, 3)

% David F. Gleich

% History
% :2011-02-16: Initial coding based on old reporting mechanism

if ~exist('topk','var') || isempty(topk), topk = min(15,length(items)); end

[ignore p] = sort(scores,'descend');
for ti=1:topk
    fprintf('%3i %8.6f %s\n',ti,scores(p(ti)), items{p(ti)});
end

