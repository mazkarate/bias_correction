function [est_moment] = compute_moments_weighted(pe_tilde, fe_tilde, weights)

    %%% Plugin estimates
    cov_mat = weighted_cov([pe_tilde, fe_tilde], weights); % Weights should sum to number of observations (NT_original)

    % Arrange to get variances and covariances in one line
    est_moment= [cov_mat(1), cov_mat(4), cov_mat(2)]; 
end