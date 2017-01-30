function [Uk,Sk,Vk,hist] = svpmc(A,k,varargin)
% SVPMC Singular Value Projection for matrix completion
% [U,S,V] = svpmc({i,j,v},k);
% [U,S,V] = svpmc({i,j,v,m,n},k);
% [U,S,V] = svpmc({N,D},k); use non-zero structure of D with valus N./D
% if D is a logical matrix, this allows zero values in N to be respected

% Copyright David F. Gleich, 2009
% University of British Columbia
% Based on svp.m by Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka
% (raghu@cs.utexas.edu)

t1=tic;

if issparse(A)
    [i,j,v]=find(A);
    [m,n]=size(A);
    i = uint(i-1);
    j = uint(j-1);
elseif iscell(A)
    if length(A)==2 % compute from A1./A2 but get non-zero structure right
        [m,n]=size(A{2});
        [i,j,v]=spratio(A{1},A{2});
        i=uint32(i-1);
        j=uint32(j-1);
    elseif length(A)==3
        i=A{1}; j=A{2}; v=A{3};
        m = max(i);
        n = max(j);
    elseif length(A)==5
        i=A{1}; j=A{2}; v=A{3};
        m=A{4}; n=A{5};
    end
end
if ~exist('m','var'), error('svpmc:invalidArgument',...
        'please provide a valid matrix input');
end
optsu=struct(varargin{:});
if isfield(optsu,'delta'), delta=optsu.delta; else delta=1/3; end
if isfield(optsu,'tol'), tol=optsu.tol; else tol=1e-3; end
if isfield(optsu,'maxit'), maxit=optsu.maxit; else maxit=500; end
if isfield(optsu,'gamma'), gamma=optsu.gamma; else gamma=0.5; end

p = numel(v)/(m*n);
eta = 1/(1+delta)/p;
nv = norm(v);

% keep everything transposed for efficency
Uk=zeros(k,m);
Sk=zeros(k,1);
Vk=zeros(k,n);
oresid = Inf;

hist = zeros(maxit,1);
fprintf('  %4s (%8s)  %8s  %8s  %8s  %8s  %8s  %8s\n', 'Iter', 'Time', ...
    'RelResid', 'DResid', 'Resid', 'RMSE', 'MAE','ResidInf');
for iter=1:maxit
    Aproj = indexprojfact(Uk,Sk,Vk,i,j);
    g = Aproj-v;
    resid = norm(g)/nv;
    resid1 = norm(g,1)/numel(v);
    residInf = norm(g,Inf);
    hist(iter)=resid;
    dt=toc(t1);
    fprintf('  %4i (%6.1f s)  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n', ...
        iter, dt, resid, abs(oresid-resid),  ... % standard residuals 
        resid*nv, resid*nv/sqrt(length(v)), resid1, residInf); % other resids
    if resid<tol || abs(oresid-resid)<tol/100, break; end
    if iter>5 && all(diff(hist(iter-4:iter))>0)
        eta = gamma*eta;
        fprintf('Residual increase detected: eta -> %g\n', eta);
    end
    oresid = resid;
    g = eta*g;
    [Uk,Sk,Vk] = projsvd(k,Uk,Sk,Vk,i,j,g);
    
end
hist = hist(1:iter);
Uk = Uk'; Vk = Vk';


function [Uk,Sk,Vk] = projsvd(k,Uk,Sk,Vk,i,j,g)
% Compute an SVD of 
% Uk*Sk*Vk' - eta*(indset(Uk*Sk*Vk',i,j) - A0)
% and return the factors, along with a residual of the new factors
% permute D and W so diagonal entries of D are sorted by proximity to sigma
boptions.tol = 1e-10 / sqrt(2);
boptions.disp = 0;
boptions.issym = 1;
m = size(Uk,2); n = size(Vk,2);
[Uk,Sk,Vk] = svdsf(@matvec,@transmatvec,[m,n],k,'L',boptions,Uk,Sk,Vk,i,j,g);
Uk = Uk';
Vk = Vk';
Sk = diag(Sk);

function y=matvec(x,Uk,Sk,Vk,i,j,g)
y = tripletmult(i,j,g,x,size(Uk,2),size(Vk,2));
y = -y;
y = y + Uk'*(Sk.*(Vk*x));

function y=transmatvec(x,Uk,Sk,Vk,i,j,g)
y = tripletmult(j,i,g,x,size(Vk,2),size(Uk,2));
y = -y;
y = y + Vk'*(Sk.*(Uk*x));
    



