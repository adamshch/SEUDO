function corrs = correlationVectorMatrix(V,M)
% compute the correlation of a vector with each column in a matrix
% V - T-length vector
% M - T x C matrix
%
% corrs - C-length matrix, each entry the correlation between the vector and a column of the matrix

% ensure correct size
assert(length(V) == size(M,1));

% ensure column
V = V(:);

% subtract means
V = V - nanmean(V);
M = bsxfun(@minus,M,nanmean(M,1));

% numerators
numers = bsxfun(@times,V,M);
numers = bsxfun(@rdivide,numers,sum(~isnan(numers))-1);

% denominators
denomV = nansum(V.^2)/(sum(~isnan(V))-1);
denomM = bsxfun(@rdivide,nansum(M.^2),sum(~isnan(M))-1);
denoms = bsxfun(@times,denomV,denomM);


% combine
corrs = nansum(numers)./sqrt(denoms);
