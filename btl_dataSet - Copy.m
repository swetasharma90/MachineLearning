function btl_dataSet(n,k,P,r)
D_set = zeros(n,n);
formatSpec = '\n%d, %d, %d ';
filename = strcat('\btl_data\dataset\Dataset_',num2str(k),'_',num2str(n));
fileID = fopen(filename,'w');
fprintf(fileID,'%d',n);
for i = 1:n
    fprintf(fileID,'\n%d',i);
end
fprintf(fileID,'\n%d',k);
if(fileID < 0)
    disp('Error');
end
for i = 1:n
    for j = (i+1):n
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
