function test_tripletmult
%%
A = [1 2 0;
     0 2 3;
     1 0 0;
     0 1 0;
     0 0 1;
     -1 0 0;
     0 -1 0;
     0 0 -1;
     1 1 1];
 A = sparse(A);
 [m,n]= size(A);
 [i,j,v] = find(A);
 ui = uint32(i-1); uj = uint32(j-1);
 ii = int32(i-1); ij = int32(j-1);
 xs = [1 0 0; 0 1 0; 0 0 1; -1 1 -1; 1 -1 1; 1.5 2.5 3.5]';
 for vi = 1:size(xs,2)
     x = xs(:,vi);
     y = A*x;
     y1 = tripletmult(i,j,v,x,m,n);
     assert(norm(y-y1)<10*eps);
     y2 = tripletmult(ui,uj,v,x,m,n);
     assert(norm(y-y2)<10*eps);
     y3 = tripletmult(ii,ij,v,x,m,n);
     assert(norm(y-y3)<10*eps);
 end
 v2 = v*exp(1);
 for vi = 1:size(xs,2)
     x = xs(:,vi);
     y = exp(1)*A*x;
     y1 = tripletmult(i,j,v2,x,m,n);
     assert(norm(y-y1)<10*eps);
     y2 = tripletmult(ui,uj,v2,x,m,n);
     assert(norm(y-y2)<10*eps);
     y3 = tripletmult(ii,ij,v2,x,m,n);
     assert(norm(y-y3)<10*eps);
 end