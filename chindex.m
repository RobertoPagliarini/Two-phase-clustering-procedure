function f = chindex(clust, X, dist)
% CHINDEX Evaluation based on the Calinski-Harabasz criterion.
%   CHINDEX(CLUST, X) computes the Calinski-Harabasz criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation/element and one
%   column per variable/feature. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Calinski-Harabasz index uses
%   the Euclidean distance between points in X.
%
%   V = CHINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Calinski-Harabasz index.
%   
%   V = CHINDEX(..., 'DISTANCE', value) computes the Calinski-Harabasz index using
%   a specified distance measure. The available built-in measures are:
%       'euc'           - Euclidean distance (the default).
%       'neuc'          - Normalized Euclidean distance.
%       'cos'           - Cosine similarity.
%       'pcorr'         - Pearson's correlation coefficient.
%
%
%
%   Reference:
%
%   T. Calinski, J. Harabasz, "A dendrite method for cluster analysis," 
%   Communications in Statistics - Theory and Methods, 
%   Vol. 3, No. 1, pp. 1-27, 1974.
%

if nargin > 2

    if strcmpi(dist,'euc')
        
        % Euclidean distance
        pfun = @eucdist;

    elseif strcmpi(dist,'neuc')
    
        % Normalized Euclidean distance
        pfun = @neucdist; 

    elseif strcmpi(dist,'cos')
    
        % Cosine similarity
        pfun = @cosdist;

    elseif strcmpi(dist,'pcorr')
    
        % Pearson's correlation coefficient
        pfun = @pcorr;
    
    else
    
        error('Unknown proximity measure');

    end

else
    
    pfun = @eucdist;

end
% ------------------------------------------------------------------------
% Validation of the clustering solution 
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<1) || K==1
    f = -inf;
    return;
end

% Intra-cluster cohesion (compactness)
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)').^2);
end
SSW = sum(sumD,1);

% Inter-cluster dispersion
Xmean = mean(X,1)';   % Global mean in X
SSB = sum(Nk .* ((feval(pfun,M',Xmean)).^2));

% CVI Computation
f = ((N-K)/(K-1))*(SSB/SSW);



function D = eucdist(A, B)
D = sqrt(bsxfun(@plus,dot(B,B,1),dot(A,A,1)')- 2*(A'*B));

function D = pcorr(A,B)
A = bsxfun(@minus,A,mean(A,2));
B = bsxfun(@minus,B,mean(B,2));
D = bsxfun(@plus,0.5*dot(B,B,1),0.5*dot(A,A,1)') - A'*B;

function D = neucdist(A, B)
A = bsxfun(@minus,A,mean(A,2));
B = bsxfun(@minus,B,mean(B,2));
A = bsxfun(@rdivide,A,std(A,0,2));
B = bsxfun(@rdivide,B,std(B,0,2));
D = sqrt(bsxfun(@plus,dot(B,B,1),dot(A,A,1)')- 2*(A'*B));

function D = cosdist(A,B)
D = bsxfun(@plus,0.5*dot(B,B,1),0.5*dot(A,A,1)') - A'*B;

