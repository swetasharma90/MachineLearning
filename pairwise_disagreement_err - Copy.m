function err= pairwise_disagreement_err(pred,P,true_ranking)
err=0;
count = 0;
for i=1:size(pred,2)
    for j=i+1:size(pred,2)
        if(find(true_ranking == pred(i))>find(true_ranking == pred(j)))
            err=err +P(pred(j),pred(i));
            count = count+1;            
        end
    end
end
end