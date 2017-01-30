function err= rank_err(pred,P)
size(pred);
P;
err=0;
for i=1:size(pred,2)
    for j=i+1:size(pred,2)
      err=err +P(pred(i),pred(j));
    end
end

end