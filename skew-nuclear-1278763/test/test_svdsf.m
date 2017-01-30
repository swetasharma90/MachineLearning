function test_svdsf
A = rand(5,3);
[U S V] = svds(A,2);
[U2 S2 V2] = svdsf(@(x) A*x, @(x) A'*x, [5 3],2);
if norm(S-S2)>eps(1)*10
    error('inaccurate svd');
end
if norm(V-V2,'fro')>eps(1)*10
    error('inaccurate svd');
end
if norm(U-U2,'fro')>eps(1)*10
    error('inaccurate svd');
end