function [est_moment] = compute_moments_group_LdM(pe_tilde, fe_tilde, ch_group, av_coworkers)

%%% Overall moment
cov_mat = fast_cov([pe_tilde, fe_tilde, av_coworkers]');

est_moment_group = zeros(numel(ch_group) - 1,6); % 6 moments: 3 variances and 3 covariances

for j= 2:numel(ch_group)
    cov_aux = fast_cov([pe_tilde(ch_group(j-1)+1:ch_group(j)), fe_tilde(ch_group(j-1)+1:ch_group(j)), av_coworkers(ch_group(j-1)+1:ch_group(j))]');
    est_moment_group(j-1,:) = [cov_aux(1), cov_aux(5), cov_aux(9),...
                    cov_aux(4), cov_aux(7), cov_aux(8) ];
end

% Arrange to get variances and covariances in one line
est_moment= [cov_mat(1), cov_mat(5), cov_mat(9),...
                    cov_mat(4), cov_mat(7), cov_mat(8) ;...
                    est_moment_group];