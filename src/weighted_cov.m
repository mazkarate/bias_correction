function Sigma2 = weighted_cov(X1, weights)

% Compute weighted mean
mu = sum(weights .* X1) / sum(weights);

% Compute weighted covariance
diff = X1 - mu;
Sigma2 = (diff' * (weights .* diff)) / (sum(weights) - 1);

