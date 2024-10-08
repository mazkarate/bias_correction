function [B_plus, B_minus] = decomp_V_LC(eps_LC, y, cluster, N_C)

ind_C1 = N_C==1; %Get indicator of cluster size equal to 1

NT = length(eps_LC(~ind_C1));

% Form auxiliary matrices
sp_eps_LC = sparse(1:NT, cluster(~ind_C1), eps_LC(~ind_C1));
sp_y      = sparse(1:NT, cluster(~ind_C1), y(~ind_C1));

% Get inner products
eps_dot_eps = full(diag(sp_eps_LC' * sp_eps_LC));
eps_dot_y   = full(diag(sp_eps_LC' * sp_y));
y_dot_y     = full(diag(sp_y' * sp_y));

% Compute the eigenvalues 
% Make sure of proper sign for extremely low values where there might be
% numerical instability. This is extremely rare
lambda_plus  = abs(0.5.*(eps_dot_y + sqrt(eps_dot_eps .* y_dot_y)));
lambda_minus = -abs(0.5.*(eps_dot_y - sqrt(eps_dot_eps .* y_dot_y)));

% Make zeros, zeros (in the extremely rare case of collinearity,
% lambda_minus = 0)
ind_zero  = abs(lambda_plus - eps_dot_eps) < 1e-12; 
lambda_minus(ind_zero) = 0;

% Compute the corresponding eigenvectors 
% For lambda_plus
a_plus = sqrt(y_dot_y ./ eps_dot_eps);
q_plus = a_plus' .* sp_eps_LC + sp_y;

% For lambda_minus
a_minus = - sqrt(y_dot_y ./ eps_dot_eps);
q_minus = a_minus' .* sp_eps_LC + sp_y;

% Normalize eigenvectors
norm_plus  = sqrt(full(diag(q_plus' * q_plus)));
norm_minus = sqrt(full(diag(q_minus' * q_minus)));

q_plus  = q_plus ./ norm_plus';
q_minus = q_minus ./ norm_minus';

% Compute B+ and B-
B_plus  = full(sum(sqrt(lambda_plus') .* q_plus, 2));
B_minus = full(sum(sqrt(abs(lambda_minus)') .* q_minus, 2));

% Now for those clusters of size 1
est_var_r = y(ind_C1).*eps_LC(ind_C1);
filt_plus = (est_var_r>=0);

est_eigenvalue = est_var_r;
est_eigenvalue(~filt_plus) = 0;
B_plus  = [B_plus; sqrt(est_eigenvalue)];

est_eigenvalue = abs(est_var_r);
est_eigenvalue(filt_plus) = 0;
B_minus  = [B_minus; sqrt(est_eigenvalue)];




% sum_test_eigenvalues = 0;
% sum_test_eigenvectors = 0;
% for i = 1:size(q_minus,2)
% 
%     % Define the matrix A
%     clust_ex = i;
%     v = eps_LC(cluster==clust_ex);
%     w = y(cluster==clust_ex);
%     A = 0.5.*(v * w' + w * v');
%     
%     % Test
%     [G , S] = eig(A);
%     
%     eig1 = max(diag(S));
%     eig2 = min(diag(S));
%     
%     %Eigenvalues
%     aux = sum(abs(lambda_plus(clust_ex) - eig1) + abs(lambda_minus(clust_ex) - eig2));
%     sum_test_eigenvalues = sum_test_eigenvalues + aux;
%     
%     d = sum(cluster==clust_ex);
%     %Eigenvectors (orthonormal vectors are not unique with respect their sign)
%     aux = sum(abs(G(:,1)) - abs(q_minus(cluster==clust_ex,clust_ex)) +...
%         abs(G(:,d)) - abs(q_plus(cluster==clust_ex,clust_ex)));
%     
%     sum_test_eigenvectors = sum_test_eigenvectors + aux;
% 
% end
% sum_test_eigenvalues
% sum_test_eigenvectors






