n=10;
A = spdiags([-ones(n,1),ones(n,1)],[-1,1],n,n);
[T,V] = slanczos_simple(A,ones(n,1),6);
norm(A*V(:,1:end-1) - V*T)
%%
% Test extension
[T3,R]=slanczos(A,ones(n,1),3);
[T5,R]=slanczos(A,ones(n,1),2,T3,R);
norm(T5-T(1:6,1:5))
%%
% Test extension with vectors
[V3,T3]=slanczos(A,ones(n,1),3,'full');
norm(A*V3(:,1:end-1)-V3*T3)
[V5,T5]=slanczos(A,ones(n,1),2,T3,V3,'full');
norm(A*V5(:,1:end-1)-V5*T5)
norm(T5-T(1:6,1:5))
norm(V5-V(:,1:6))
%%
[V6,T6]=slanczos(A,ones(n,1),3,T3,V3,'full');
norm(A*V6(:,1:end-1)-V6*T6)