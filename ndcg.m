function undcg =  ndcg(true_score,pred_score,pred_rank)
pred_score=sort(pred_score,'descend');
dcg = 0.0; 
rank = 1;
for i = 1:size(pred_score,2)
    dcg = dcg+computeDCG(pred_score(1,i), pred_rank(1,i)); 
    %rank = rank+1;
end
% score_vec = true_score;
% idcg = 0.0; 
% rank = 1;
% score_vec = sort(score_vec,'descend');
% for i = 1:size(score_vec,2)
%     idcg = idcg+computeDCG(score_vec(1,i), rank); 
%     rank = rank+1;
% end
% dcg
% idcg
%  undcg = dcg / idcg;
undcg=dcg;
    function dcg= computeDCG(rel,rank) 
        dcg = ((2^rel) - 1.0) /( log(rank + 1) /log(2)); 
    end
        
end