function [U S V]=svp(omega,b,k,m,n)

%get indices back for making matrix from vector of non-zero values
[I, J] = ind2sub([m,n], omega);

t = 1;
X = zeros(m,n);

t_norm = norm(X(omega) - b);

epsilon = 0.0001;
max_iter = 100
iter = 0;
while(t_norm > epsilon & t<max_iter)
    k_svd=X(omega) - (1/t)*(X(omega) -b);
    
    %making matrix to perform svd
    svd_X=zeros(m,n);
    for i=1:size(b)
        svd_X(I(i),J(i))=k_svd(i);
    end
    
    [U,S,V] = svd(svd_X);
    %k-rank approximation
    X = U(1:m,1:k)*S(1:k,1:k)*V(1:n,1:k)';
    %X = U*S*V';
    t = t+1;
    t_norm = norm(X(omega) - b);
    iter = iter+1;
end
end
