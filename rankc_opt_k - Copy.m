function rankc_opt_k()
min_err_index = zeros(1,20);
n_log_base_e = zeros(1,20);
%for n = 5:20
n=50;
    i=1;
    k_start =1
    k_end=50
    [P,true_r,weight_vec] = generate(n,k_start,k_end);
    true_err=rank_err(true_r,P)
    for k=1:50 % if u change the rate of change of 'K' then please change in min_err_index also
%            if(k==122)
%                pd_err(i)=0;
%                i=i+1;
%            continue;
%           end
k
        file = strcat('btl_data\dataset\Dataset_',num2str(k),'_',num2str(n));
        [index score]=nnm(file);
        %index=svm_aggregation(file);
        file
       % [index score]=rankcentrality(file);
       % err(i)=rank_err(index,P);
       % pd_err(i) = pairwise_disagreement_err(index,P,true_r)
    index
    score
       pd_err(i) = ndcg(weight_vec,score,index)
        i=i+1;
    end   
%     [M,I] = min(pd_err)
%     min_err_index(n) = 10+(I-1);
%     n_log_base_e(n) = n*log(n);
%end

%subplot(1,2,1)

plot(1:50,pd_err)
% xinitial=1;
% xfinal=50;
% yinitial=0;
% yfinal=0.2;
% axis([xinitial xfinal yinitial yfinal]);
%,1:20,n_log_base_e)
title('n=50 and k = 50 to 100 vs error');
% legend('optimal_k','nlogn')
% title('n vs optimal k');
% subplot(1,2,2)
% plot(1:20,min_err_index,1:20,n_log_base_e)
% legend('optimal_k','nlogn')
% title('n vs optimal k');


end
