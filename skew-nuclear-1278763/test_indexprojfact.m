function test_indexprojfact
%%
A = [0 1 0
     1 0 0
     1 1 1];
[i,j] = find(A);
ui = uint32(i-1); uj = uint32(j-1);
ii = int32(i-1); ij = int32(j-1);

randnseed = randn('seed');
randn('seed',1);
k = 2;
U = rand(size(A,1),k);
S = rand(k,1);
V = rand(size(A,2),k);

v = nonzeros(A.*(U*diag(S)*V'));
vt = indexprojfact(U',S,V',i,j);
assert(norm(v-vt)<10*eps);
vt = indexprojfact(U',S,V',ii,ij);
assert(norm(v-vt)<10*eps);
vt = indexprojfact(U',S,V',ui,uj);
assert(norm(v-vt)<10*eps);

