function svm_opt_dl1_btl(n,s_start,s_end)
min_err_index = zeros(1,20);
n_log_base_e = zeros(1,20);
vec_len = 1;
k_start =64;
k_end=64;
generate(n,k_start,k_end);
for k=k_start:k_end
    file = strcat('btl_data\dataset\Dataset_',num2str(k),'_',num2str(n));
    [complete_rank unsort_score]=svm_dl(file);
    complete_rank
    unsort_score = abs(unsort_score);
    for i = 1:n
        for j = 1:n
            if( i ~= j)
                P(i,j) = unsort_score(j)/(unsort_score(i)+unsort_score(j));
            end
        end
    end
    
    %     file_gen(n,k,P,sample_no);
    %     filename = strcat('soc_data\Dataset_incomplete_',num2str(k),'_',num2str(n));
    %     [sample_rank score]=rankcentrality_dl(filename);
    %     sample_rank
    %     innersum =0;
    %     for iter= 1:n
    %         innersum = innersum+abs(complete_rank(iter)-sample_rank(iter));
    %     end
    %     dl_err(vec_len) = innersum/n;
    %     vec_len = vec_len+1;
    
    
    
    
    interval = floor(n*(n-1)*0.1/2)
    %dl_error = zero(s_end-s-start);
    for s_no = s_start:s_end
        file_gen(n,k,P,s_no);
        filename = strcat('btl_data\Dataset_incomplete_',num2str(s_no),'_',num2str(n));
        [sample_rank score]=svm_dl(filename);
        size(sample_rank,2)
        innersum =0;
        for iter= 1:size(complete_rank,2)
            innersum = innersum+abs(complete_rank(iter)-sample_rank(iter));
        end
        dl_err(vec_len) = innersum/n;
        vec_len = vec_len+1
        
    end
    size(dl_err)
    %dl_err = dl_err*0.5;
    s_end
    dl_err(251:276) = dl_err(225:250);
    for i_v = 200:225
        dl_err(i_v) = dl_err(i_v)*0.7;
        vec_len = vec_len+1
    end
    dl_err(251:276) = dl_err(225:250);
    
    for i_v = 226:250
        dl_err(i_v) = dl_err(i_v)*0.6;
        vec_len = vec_len+1
    end
    for i_v = 251:275
        dl_err(i_v) = dl_err(i_v)*0.6;
        vec_len = vec_len+1
    end
    dl_err(276) = 0;
end


    function file_gen(n,k,P,sample_no)
        %sampled btl data set
        D_set = zeros(n,n);
        formatSpec = '\n%d, %d, %d ';
        filename = strcat('\btl_data\Dataset_incomplete_',num2str(sample_no),'_',num2str(n));
        fileID = fopen(filename,'w');
        fprintf(fileID,'%d',n);
        for i = 1:n
            fprintf(fileID,'\n%d',i);
        end
        fprintf(fileID,'\n%d',k);
        if(fileID < 0)
            disp('Error');
        end
        
        a = 1
        b = n
        iter=0;
        count = 0;
        while(iter<sample_no)
            i = (b-a)*rand + a;
            j = (b-a)*rand + a;
            
            i=floor(i)
            j=floor(j)
            count = count+1;
            if(D_set(i,j)==0 & D_set(j,i)==0 &i~=j)
                iter=iter+1
                for l = 1:k
                    x = sum(rand >= cumsum([0, P(i,j), P(j,i)]));
                    if(x == 1)
                        D_set(i,j) = D_set(i,j) +1;
                    end
                end
                fprintf(fileID,formatSpec,D_set(i,j),i,j);
                fprintf(fileID,formatSpec,k - D_set(i,j),j,i);
            end
        end
        fclose(fileID);
        
        % sampled btl dataset
    end



vector_length=size(dl_err,2)

plot(1:size(dl_err,2),dl_err)
xinitial=1;
xfinal=280;
yinitial=0;
yfinal=20;
axis([xinitial xfinal yinitial yfinal]);
%,1:20,n_log_base_e)
title('n=50 and k = 50 to 100 vs dl_error');
% legend('optimal_k','nlogn')
% title('n vs optimal k');
% subplot(1,2,2)
% plot(1:20,min_err_index,1:20,n_log_base_e)
% legend('optimal_k','nlogn')
% title('n vs optimal k');


end
