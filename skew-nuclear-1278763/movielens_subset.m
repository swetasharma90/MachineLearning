%% Build a subset of movielens
A=readSMAT('/home/dgleich/data-new/movielens/100k/movielens.smat');
movies = readList('/home/dgleich/data-new/movielens/100k/movielens.names');

% the goal is to find 20 users and 8 movies such that all users have rated
% at least 5 things, and all movies have at least 5 ratings.
rand('seed',0); % ensure repeatability
nu = 20;
nm = 8;
k = 5;
mr = 3;
for i=1:200
    us = 1:(nu+4*i);
    ms = 1:(nm+2*i);
    As = A(us,ms);
    A1 = spones(As);
    numvalidusers = sum(sum(A1,2)>=k);
    numvalidmovies = sum(sum(A1,1)>=k);
    if numvalidusers>=nu && numvalidmovies>=nm,
        us = us(sum(A1,2)>=k);
        ms = ms(sum(A1,1)>=k);
        if length(us)>nu,
            p = randperm(length(us));
            us = us(p(1:nu));
        end
        if length(ms)>nm,
            p = randperm(length(ms));
            ms = ms(p(1:nm));
        end
        
        break;
    end
end
% add Taxi Driver in
ms = [ms 23];
curuser = nu+4*i+1; % make sure it isn't someone we saw
while 1
    As = A(us,ms);
    A1 = spones(As);
    movies_small = movies(ms);
    invalidusers = sum(A1,2)<mr;
    if sum(invalidusers) == 0
        break
    end
    us(invalidusers) =  [];
    us = [us curuser:curuser+sum(invalidusers)-1];
    curuser = curuser + sum(invalidusers);
end


%