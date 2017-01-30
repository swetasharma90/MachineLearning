%% Reproduce figure 4 from Meka, Jain and Dhillon
% 
ns = 10:500:5000;
k = 2;
p = 0.1;
dts = zeros(size(ns));
for ni=1:length(ns)
    n = ns(ni);
    As = sprand(n,n,p);
    A = randn(n,2)*randn(2,n);
    [i,j,v] = find((As>0).*A);
    tic
    [Uk,Sk,Vk,t,err]=svp([i,j],v,n,n,k);
    dts(ni) = toc;
end
