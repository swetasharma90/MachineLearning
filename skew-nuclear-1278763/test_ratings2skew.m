function test_ratings2skew
%% single user tests
A=[1 0 4];
[Y,N] = ratings2skew(sparse(A),'am');
Ytrue = [0 0 -3; 0 0 0; 3 0 0];
Ntrue = [0 0 1; 0 0 0; 1 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 0 4];
[Y,N] = ratings2skew(sparse(A'),'am',[],[],2);
Ytrue = [0 0 -3; 0 0 0; 3 0 0];
Ntrue = [0 0 1; 0 0 0; 1 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 0 4];
[Y,N] = ratings2skew(sparse(A),'gm');
Ytrue = [0 0 log(1)-log(4); 0 0 0; log(4)-log(1) 0 0];
Ntrue = [0 0 1; 0 0 0; 1 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 0 4];
shift = exp(1);
[Y,N] = ratings2skew(sparse(A),'gm',[],shift);
Ytrue = [0 0 log(1+shift)-log(4+shift); 0 0 0; log(4+shift)-log(1+shift) 0 0];
Ntrue = [0 0 1; 0 0 0; 1 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 0 4];
[Y,N] = ratings2skew(sparse(A),'bc');
Ytrue = [0 0 -1; 0 0 0; 1 0 0];
Ntrue = [0 0 1; 0 0 0; 1 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 0 4];
[Y,N] = ratings2skew(sparse(A),'bc',[],[],[],'equiv',5);
Ytrue = [0 0 0; 0 0 0; 0 0 0];
Ntrue = [0 0 1; 0 0 0; 1 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 0 4];
[Y,N] = ratings2skew(sparse(A),'bc',[],[],[],'equiv',1);
Ytrue = [0 0 -1; 0 0 0; 1 0 0];
Ntrue = [0 0 1; 0 0 0; 1 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))


A=[1 1 4];
[Y,N] = ratings2skew(sparse(A),'sb',[],[],[]);
Ytrue = [0 0 -1; 0 0 -1; 1 1 0];
Ntrue = [0 0 1; 0 0 1; 1 1 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 1 4];
[Y,N] = ratings2skew(sparse(A),'sb',[],[],[],'equiv',5);
Ytrue = [0 0 0; 0 0 0; 0 0 0];
Ntrue = [0 0 0; 0 0 0; 0 0 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 1 4];
[Y,N] = ratings2skew(sparse(A),'sb',[],[],[],'equiv',2);
Ytrue = [0 0 -1; 0 0 -1; 1 1 0];
Ntrue = [0 0 1; 0 0 1; 1 1 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A = [1 2; 2 1; 1 1];
[Y,N] = ratings2skew(A,'lo',[],[],[],'pseudo',0);
Ytrue = [0 0; 0 0]; Ntrue = [0 3; 3 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A = [1 2; 2 1; 1 1];
[Y,N] = ratings2skew(A,'lo',[],[],[],'pseudo',1);
Ytrue = [0 0; 0 0]; Ntrue = [0 3; 3 0];
assert(norm(Ytrue-Y)<10*eps && isequal(N,Ntrue))

A=[1 0 4];
[Y,N] = ratings2skew(A,'lo',[],[],[],'pseudo',1);
Ytrue = [0 0 log(2/1); 0 0 0; log(1/2) 0 0];
Ntrue = [0 0 1; 0 0 0; 2 0 0];

%% Conversion to sparse
A=[1 1 4];
[Y1,N1] = ratings2skew(A,'sb',[],[],[],'equiv',2);
[Y,N] = ratings2skew(sparse(A),'sb',[],[],[],'equiv',2);
assert(norm(Y1-Y)<10*eps && isequal(N,N1))
    
%%
% test ratings to skew
% A=readSMAT('/home/dgleich/data-new/movielens/100k/movielens.smat');
% [Y,N]=ratings2skew_mex(A,At,4);
% [Y1,N1]=ratings2skew(A);
% % find everyone who rated movies 1 and 2
% m1_users=find(A(:,1));
% m2_users=find(A(:,2));
% both_users = intersect(m1_users,m2_users);
% rats = full([A(both_users,1) A(both_users,2)]);
% sum(rats(:,1)>rats(:,2))
