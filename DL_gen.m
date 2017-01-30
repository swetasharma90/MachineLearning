function DL_gen(file,n,sample_no,k_start,k_end)

     
fid = fopen(file, 'r');
tline = fgetl(fid);
num = str2num(tline);
ctr = 0;
soc=zeros(num,num);
while(~feof(fid))
    str = fgetl(fid);
    ctr = ctr + 1;
    if(ctr>(num+1))
        index = ctr-num-1;
        B(index,:) = regexp(str, '\,', 'split');
        if(str2double(B(index,1))~=0)
        soc_mat(index,:)=str2double(B(index,:));       
        end
    end
end
r = zeros(num,num);
for i= 1:size(soc_mat,1)
    for j= 2:size(soc_mat,2)
        for l= (j+1):size(soc_mat,2)
            r(soc_mat(i,l),soc_mat(i,j)) = r(soc_mat(i,l),soc_mat(i,j)) + soc_mat(i,1);
        end
    end    
end
for i = 1:num
    for j = 1:num
        if( i ~= j)
            P(i,j) = r(i,j)/(r(i,j)+r(j,i));
        end
    end
end
     
P
for k=k_start:k_end
    file_gen(n,k,P,sample_no);
end


    function file_gen(n,k,P,sample_no)
        %sampled btl data set
        D_set = zeros(n,n);
        formatSpec = '\n%d, %d, %d ';
        filename = strcat('\btl_data\dataset\Dataset_incomplete_',num2str(k),'_',num2str(n));
        fileID = fopen(filename,'w');
        fprintf(fileID,'%d',n);
        for i = 1:n
            fprintf(fileID,'\n%d',i);
        end
        fprintf(fileID,'\n%d',k);
        if(fileID < 0)
            disp('Error');
        end
        
        a = 1;
        b = n;
        iter=0;
        while(iter<sample_no)
            i = (b-a)*rand + a;
            j = (b-a)*rand + a;
            i
            j
            i=floor(i);
            j=floor(j);
            if(D_set(i,j)==0 & D_set(j,i)==0)
                iter=iter+1;
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
end