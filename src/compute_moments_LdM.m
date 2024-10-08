function [est_moment] = compute_moments_LdM(pe_tilde, fe_tilde, av_coworkers)

    %%% Overall moment
    cov_mat = fast_cov([pe_tilde, fe_tilde, av_coworkers]');

    % Arrange to get variances and covariances in one line
    est_moment= [cov_mat(1), cov_mat(5), cov_mat(9),...
                    cov_mat(4), cov_mat(7), cov_mat(8) ];
end