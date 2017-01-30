function [Y,N]=ratings2skew(A,type,fullflag,shift,dim,varargin)
% RATINGS2SKEW Convert ratings data to a skew-symmetric matrix
%
% Y = ratings2skew_simple(A,type,fullflag,shift,dim,...)
% A : a ratings matrix which is usually nusers x nitems.
% type : a string to pick the agglomeration type of the ratings
%   each option has a short or long form.
%   'am' - arithmetic mean of score differences ('arithmean')
%   'gm' - geometric mean of score rations ('geomean')
%   'bc' - binary comparison ('binary')
%   'sb' - strict binary ('strictbinary') (equivalent to bc without ties)
%   'lo' - log odds ('logodds')
% full : a binary flag indicating the result should be fairly dense
%   and so we should treat it as a full matrix
%   0 (default) - the matrix should be sparse
%   1 - the matrix should be full
% shift : a positive scalar shift to accomodate zero ratings
% dim : a variable for the user dimension of the ratings matrix
%   1 (default) - each user is a row of the matrix
%   2  - each user is a column of the matrix.  
% At the moment, this function only does full computations.
% You can also specify a few extra options:
% 'equiv' - specifies an equivalence region for ratings.  In other
%   words, it treats differences smaller than equiv as 0.  This is used
%   used for the two 
% 'pseudo' - a count to add pseudocount smoothing to the log-odds ratio
%   to correct for 0 or 1 probabilities (which cause issues in the 
%   numerator or denominator of the log odds function.  (default = 1, 
%   which implies that we assume each baseline probability is 1/2.
%   A good choice for this is the minimum number of ratings on any pair
%   of items.

% David F. Gleich, Copyright 2009
% University of British Columbia

% History
% :2009-11-21: Initial coding
% :2009-12-05: Added log-odds

if ~exist('fullflag','var') || isempty(fullflag), fullflag=0; end
if ~exist('shift','var') || isempty(shift), shift=0; end
if ~exist('dim','var') || isempty(dim), dim=1; end
optsu = struct(varargin{:});
if isfield(optsu,'equiv'), equiv=optsu.equiv; else equiv=0; end
if isfield(optsu,'pseudo'), pseudo=optsu.pseudo; else pseudo=1; end

% check if A is sparse
if ~issparse(A), A = sparse(A); end

if dim==1, A = A'; end; % TODO special case this transpose in the future
% Now users correspond with columns

ntype = [];
switch type
    case {'am','arithmean'}, ntype=1;
    case {'gm','geoman'}, ntype=2;
    case {'bc','binary'}, ntype=3;
    case {'sb','strictbinary'}, ntype=4;
    case {'lo','logodds'}, ntype=5; 
        if pseudo==0, warning('ratings2skew:mayBeInvalidRatio',...
                'the ratio Y./N may not be valid without pseudo-counts');
        end
    otherwise error('ratings2skew:invalidType','type %s is not valid',type);
end
[Y,N]=ratings2skew_mex(A,ntype,shift,equiv,pseudo);
switch type
    case {'lo','logodds'}
        % N here counts the number of comparisons.  
        % Y(i,j) counts the number of times i beat j
        % so the return is Y = log(Y./Y').*N over the non-zeros of N so that
        % Y./N gives the correct value.
        % Y = splogratio3(Y,Y',N); 
        % THIS IS NOW IMPLEMENTED INSIDE THE MEX FILES
end

% function Y=splogratio3(Y1,Y2,N)
%     % compute (Y1./Y2).*N over the non-zeros of N, which should
%     % be a super-set of the non-zeros of Y1, Y2 (if not, it just
%     % drops those elements.
%     m = size(N,1); n = size(N,2);
%     [i,j,nval]=find(N); % get the full non-zero structure
% 
%     % build a map between (i,j) -> non-zero element in v2
%     Ni=sparse(i,j,1:length(nval),m,n);
%     y1i = nonzeros((Y1~=0).*Ni); % find the indices of non-zeros in A1
%     y1v = nonzeros(Y1.*spones(N));
%     y1 = zeros(size(nval));
%     y1(y1i) = y1v;
%     y2i = nonzeros((Y2~=0).*Ni); % find the indices of non-zeros in A1
%     y2v = nonzeros(Y2.*spones(N));
%     y2 = zeros(size(nval));
%     y2(y2i) = y2v;
%     
%     
%     Y = sparse(i,j,log(y1./y2).*nval,m,n);
%     
%     
