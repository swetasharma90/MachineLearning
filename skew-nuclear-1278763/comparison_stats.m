function S=comparison_stats(N)
S.nnz = nnz(N);
S.num_movies=sum(sum(N)>0);
S.num_missing_movies=size(N,1)-S.num_movies;
S.p = nnz(N)/numel(N);
S.max_comparisons = max(max(N));
S.min_comparisons = min(nonzeros(N));
