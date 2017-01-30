%% 

A1 = sparse(6,4);
A1(1,2) = 2;
A1(2,3) = 4;
A1(4,4) = 3;
A1(4,1) = -2;
A2 = sparse(6,4);
A2(1,2) = 5;
A2(2,1) = 4;
A2(2,3) = 5;
A2(4,4) = 6;
A2(4,1) = 2;
A2(3,3) = 2;
[Uk,Sk,Vk,hist] = svpmc({A1,A2},2);

%%
addpath('svp-1.0');
A = sprand(5,3,0.1*5);
[i,j,v] = find(A);
[m,n] = size(A);
[Uk,Sk,Vk,hist]=svpmc({i,j,v,m,n},2); semilogy(hist)
[U2,S2,V2,t,err]=svp([i,j],v,m,n,2);
[norm(Uk-U2) norm(Sk-S2) norm(Vk-V2)]

