function dataSet_soc(file,k)
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
D_set = zeros(num,num);
formatSpec = '\n%d, %d, %d ';
filename = strcat('soc_data\Dataset_from_Soc',num2str(k),'_',num2str(num))
fileID = fopen(filename,'w');
fprintf(fileID,'%d',num);
for i = 1:num
    fprintf(fileID,'\n%d',i);
end
fprintf(fileID,'\n%d',k);
for i = 1:num
    for j = (i+1):num
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
end
