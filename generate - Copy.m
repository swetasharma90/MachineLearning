function [P,index,r_sort]=generate(n,k_start,k_end)

r = rand(1,n);
for i = 1:n
    for j = 1:n
        if( i ~= j)
            P(i,j) = r(j)/(r(i)+r(j));
        end
    end
end

for k=k_start:k_end
btl_dataSet(n,k,P,r);
end
[r_sort index] = sort(r,'descend');
filename = strcat('\btl_data\truerank\True_ranking_','_',num2str(n));
dlmwrite(filename,index);
dlmwrite(filename, r_sort, '-append');


end
