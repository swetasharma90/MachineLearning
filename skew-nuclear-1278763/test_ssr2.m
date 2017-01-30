%% Setup
normdiffshift = @(x,y) norm((x(:) - x(1)) - (y(:) - y(1)));
types = {'am','gm','bc','sb','lo'};

%% Seeded example
n = 3;
s = 1:n;
%Y - (ge^T - eg^T)
Y =  s(:)*ones(1,n) - ones(n,1)*s(:)';

e = ones(n,1);
n = 3;
A = kron(e,speye(n)) - kron(speye(n),e);
g = A\Y(:);
assert(normdiffshift(s,g)<100*eps, 'failed small seeded example');


% test complete flag, but it does't matter here
[AtA,b] = ssr2(Y,1);
g2 = pinvsolve(AtA,b);
assert(normdiffshift(s,g2)<100*eps, 'failed small seeded example');

AtA2 = AtA(1:end-1,1:end-1); b2=b(1:end-1);
g3 = AtA2\b2; g3(end+1) = 0;
assert(normdiffshift(s,g3)<100*eps, 'failed small seeded example');

% test sparse input
[AtA,b] = ssr2(sparse(Y));
g2 = pinvsolve(full(AtA),b);
assert(normdiffshift(s,g2)<100*eps, 'failed small seeded example');

[i,j,v] = find(Y);
[AtA,b] = ssr2({i,j,v,size(Y,1), size(Y,2)});
g2 = pinvsolve(full(AtA),b);
assert(normdiffshift(s,g2)<100*eps, 'failed small seeded example');



%% Try a small real example
A = [
     5     3     0     0     0     4     0     0
     4     5     5     3     0     1     1     1
     0     5     4     5     4     0     0     1
     0     4     5     0     5     0     4     1
     5     3     0     0     0     0     1     1
     0     5     0     4     4     0     0     0
     0     2     0     4     0     3     1     0
     3     5     0     5     0     4     0     0
     0     5     4     5     4     0     0     1
     4     0     5     0     3     0     0     0
     0     4     5     0     5     0     4     1
     0     5     0     2     4     1     0     0
     0     5     4     0     3     0     0     0
     4     3     0     3     2     0     0     0
     0     5     0     4     2     2     0     0
     0     3     5     0     3     0     0     0
     0     5     0     4     3     0     0     0
     0     4     0     3     0     2     0     0
     4     3     4     0     0     2     0     0
     0     3     0     2     3     1     0     0
];
movies = {
    'Dr. Strangelove or How I Stopped  (1964)'
    'Groundhog Day (1993)'
    'As Good As It Gets (1997)'
    'Home Alone (1990)'
    'Grumpier Old Men (1995)'
    'Brady Bunch Movie, The (1995)'
    'Leave It to Beaver (1997)'
    'Mr. Magoo (1997)'};


for ti=1:length(types)
    [Y,N] = ratings2skew(A,types{ti});

    Yr = full(spratio(Y,N));

    n = size(Y,1);
    e = ones(n,1);
    Asys = kron(e,speye(n)) - kron(speye(n),e);

    g = Asys\Yr(:);
    [ignore p] = sort(g,'descend');
    movies(p)

    Z = ones(n,1)*g(:)' -g(:)*ones(1,n);

    [AtA,b] = ssr2(spratio(Y,N),1);
    g2 = pinvsolve(full(AtA),b);
    assert(normdiffshift(g,g2)<100*eps, 'failed small example');
    
    % test sparse input
    if nnz(N)==numel(N)-size(N,1)
        [i,j,v] = spratio(Y,N);
        [AtA,b] = ssr2({i,j,v,size(Y,1), size(Y,2)});
        g2 = pinvsolve(full(AtA),b);
        assert(normdiffshift(g,g2)<100*eps, 'failed small example');
    end
end

%% Planted example with noise
randn('seed',0);
n = 5;
s = 1:n; s = s(:);
sigma2 = 1/sqrt(2);
nusers = 25;
A = repmat(s',nusers,1);
A = A + randn(size(A))*sqrt(sigma2);
A(A>5) = 5;
A(A<1) = 1;
A = round(A);

for ti=1:length(types)
    [Y,N] = ratings2skew(A,types{ti});
    Yr = full(spratio(Y,N));

    n = size(Y,1);
    e = ones(n,1);
    Asys = kron(e,speye(n)) - kron(speye(n),e);

    g = Asys\Yr(:);
    Z = ones(n,1)*g(:)' -g(:)*ones(1,n);

    [AtA,b] = ssr2(spratio(Y,N),1);
    g2 = pinvsolve(full(AtA),b);
    assert(normdiffshift(g,g2)<100*eps, 'failed noise example');
    
    if nnz(N)==numel(N)-size(N,1)
        [i,j,v] = spratio(Y,N);
        [AtA,b] = ssr2({i,j,v,size(Y,1), size(Y,2)});
        g2 = pinvsolve(full(AtA),b);
        assert(normdiffshift(g,g2)<100*eps, 'failed small example');
    end
end

