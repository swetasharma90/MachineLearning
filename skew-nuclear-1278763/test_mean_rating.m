function test_mean_rating

R = [1 1
     2 0];
mu = mean_rating(R);
assert(mu(1) == 1.5);
assert(mu(2) == 1);
mu = mean_rating(R,2);
assert(mu(1) == 1);
assert(mu(2) == 2);
R = zeros(10,4);
mu = mean_rating(R);
mu
assert(all(mu==0));
assert(length(mu)==4);