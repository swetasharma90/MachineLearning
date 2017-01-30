function nnm(file)
rank = 2; %rank of the optimal matrix
fid = fopen(file, 'r');
tline = fgetl(fid);

num = str2num(tline);
ctr = 0;

%Reading File
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
        k = str2num(str)
    end    
end

%constructing Y
for i=1:size(a,1)
    for j=1:size(a,2)
        if i~=j & a(j,i) ~=0 & a(i,j) ~=0
            Y(i,j) = log(a(i,j)/a(j,i));
        end
        if a(j,i) ==0
            Y(i,j) =0;
        end
        if a(i,j) ==0
            Y(i,j) =0;
        end
    end
end

%dummy File
%Y=[0 0.2 0.3 0.4;-0.2 0 0.6 0.7;-0.3 -0.6 0 0.8;-0.4 -0.7 -0.8 0];

Omega = find(Y);
m = size(Y,1);
n = size(Y,2);

[U,S,V] = svp(Omega, Y(Omega), rank,m,n); %Run SVP

s=U*S*V'*ones(n,1)/n

end

