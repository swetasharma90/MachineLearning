function [A0,ss,cn] = coresplit(A,k)
% CORESPLIT Split a matrix into a smaller matrix based on core numbers
%
% A0=coresplit(A,k) returns all vertices of A or [0 A; A' 0] with core
% number larger than or equal to k.

if size(A,1)~=size(A,2)
    B = spaugment(A,0);
    cn = corenums(B);
    f = cn>=k;
    i1 = f(1:size(A,1));
    i2 = f(size(A,1)+1:end);
    A0 = A(i1,i2);
    ss = {i1,i2};
else
    cn = corenums(A);
    f = cn>=k;
    ss = f;
    A0 = A(ss,ss);
end