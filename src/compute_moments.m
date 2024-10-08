function [est_moment] = compute_moments(pe_tilde, fe_tilde)

    %%% Overall moment
    cov_mat = fast_cov([pe_tilde, fe_tilde]');

    % Arrange to get variances and covariances in one line
    est_moment= [cov_mat(1), cov_mat(4), cov_mat(2)];
end