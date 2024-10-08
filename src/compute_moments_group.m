function [est_moment] = compute_moments_group(pe_tilde, fe_tilde, ch_group)

%%% Overall moment
cov_mat = fast_cov([pe_tilde,fe_tilde]');

est_moment_group = zeros(numel(ch_group) - 1,3);

for j= 2:numel(ch_group)
    cov_aux = fast_cov([pe_tilde(ch_group(j-1)+1:ch_group(j)), fe_tilde(ch_group(j-1)+1:ch_group(j))]');
    est_moment_group(j-1,:) = [cov_aux(1), cov_aux(4), cov_aux(2)];
end

% Arrange to get variances and covariances in one line
est_moment= [cov_mat(1), cov_mat(4), cov_mat(2); est_moment_group];