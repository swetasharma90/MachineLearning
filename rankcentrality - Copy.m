function [index, p_sort]=rankcentrality(file)

fid = fopen(file, 'r');
tline = fgetl(fid);

num = str2num(tline);
a = zeros(num,num);
ctr = 0;

while(~feof(fid))
    str = fgetl(fid);
    ctr = ctr + 1;
    if(ctr>(num+1))
        index = ctr-num-1;
        B(index,:) = regexp(str, '\,', 'split');
        win=str2double(B(index,1));
        i=str2double(B(index,2));
        j=str2double(B(index,3));
        a(i,j)=win/k;
    end
    if(ctr == (num+1))
       k = str2num(str); 
    end
    
end
d=sum(a~=0,2);
d_max=max(d);
P=a/d_max;
for i=1:num
    P(i,i)=1-sum(P(i,:));
end
[V,D] = eig(P);
[maxval,index] = max(diag(D));
maxVector = V(:,index);
p=(1/num)*ones(1,num);

 while(true)
     p_t=p*P;
    if(p_t==p)
        break;
    end
    p=p_t;
end
[p_sort,index] = sort(p_t,'descend');
index;
fclose(fid);
end