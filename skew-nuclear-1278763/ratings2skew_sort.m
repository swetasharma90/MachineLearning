function [S,Sc] = ratings2skew(A)
% RATINGS2SKEW Convert from ratings data into skew symmetric data
% S = ratings2skew(A) returns a skew-symmetric matrix from a set
% of ratings where any movie -- say mi -- with a rating below another 
% movie -- say mj -- from the same user causes an entry 
% sparse([mi;mj],[mj;mi],[-1;1]) to be added to the matrix S
% [S,Sc] = ...
% also returns the counts of the elements, this allows normalizing
% by the total number of evaluated cells
% In the matrix A, users are rows, movies are columns

nmovies=size(A,2);
nusers=size(A,1);
bufsize=2e7;
S = sparse(nmovies,nmovies);
Sc = S;
At = A'; % get access to users
pi = zeros(bufsize,1); pj=zeros(bufsize,1);
pos = 1;
for i=1:nusers
  [mi ignore ri] = find(At(:,i)); % user i has ratings ri on movies mI
  [ignore p] = sort(ri); % sort the ratings (optimize with bucket-sort)
  mi = mi(p); ri=ri(p);
  nrat = length(mi);
  if nrat^2 + pos > bufsize
    Sd=update_S(nmovies,pi,pj,pos);
    S = S - Sd + Sd'; Sc = Sc+Sd+Sd'; 
    pos=1;
  end
  % given nrat entries, we are going to put in at most nrat^2 
  % skew-symmetric entries
  currat=ri(1);
  lastrat=1;
  for k=2:nrat
    if ri(k)~=currat
      % then movies 1:k-1 are WORSE than all the movies that come after
      for mi1=lastrat:k-1
        for mi2=k:nrat
          pi(pos)=mi(mi1); pj(pos)=mi(mi2); pos=pos+1;
        end
      end
      currat=ri(k);
      lastrat=k;
    end
  end
end
if pos>0
  Sd = update_S(nmovies,pi,pj,pos);
  S = S - Sd + Sd'; Sc = Sc+Sd+Sd'; 
  pos=1;
end

function Sd=update_S(n,pi,pj,pos)
% pi is the movie that is worse than pj
Sd = sparse(pi(1:pos-1),pj(1:pos-1),1,n,n);
