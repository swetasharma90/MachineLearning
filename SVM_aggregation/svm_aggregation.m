function sigma=svm_aggregation(file)

fid = fopen(file, 'r');
tline = fgetl(fid);

num = str2num(tline);
ctr = 0;
%reading file
while(~feof(fid))
    str = fgetl(fid);
    ctr = ctr + 1;
    if(ctr>(num+1))
        index = ctr-num-1;
        B(index,:) = regexp(str, '\,', 'split');
        win=str2double(B(index,1));
        i=str2double(B(index,2));
        j=str2double(B(index,3));
        P_hat(i,j)=win/k;
    end
    if(ctr == (num+1))
        k = str2num(str);
    end
    
end
% generating S_p dataset
c=1;
for i=1:size(P_hat,2)
    for j=1:size(P_hat,2)
        if(i<j)
            S_p(c,:)= [(P_hat(:,i)-P_hat(:,j))' sign(P_hat(j,i)-P_hat(i,j))];
            if(S_p(c,end)==0)
                S_p(c,end)=1;
            end
            c=c+1;
        end
    end
end
n=size(P_hat,1);

%linearly seperable checking
train_data=S_p(:,1:end-1);
train_label=(-1).*S_p(:,end);

for i=1:size(train_data,1)
    A(i,:)=train_data(i,:).*train_label(i,:);
end
Q=eye(size(S_p,2)-1);
a=zeros(size(S_p,1),1);
q=zeros(size(S_p,2)-1,1);
w0=zeros(size(S_p,2)-1,1);
options = optimset('Algorithm','interior-point-convex','display','off','MaxIter',10000);
[w, objective,exitflag] = quadprog(Q,q,A,a,[],[],[],[],[],options);
C=1/size(S_p,1);
if(exitflag>0)
    %hard margin
    model=SVM_learner('h',train_data, train_label, C, 'linear', 1);
else
    %soft margin
    model=SVM_learner('s',train_data, train_label, C, 'linear', 1);
end
model.w
for i=1:n
    f(i)=0;
    for k=1:n
        f(i)=f(i)+model.w(k)*P_hat(k,i);
    end
end
[f_sort,sigma] = sort(f,'descend');
sigma %argsort f

end