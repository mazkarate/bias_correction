function [P_hat, M_hat] = est_M_P(n_lev, NT_mover, X, xx, L, cluster, mytol)

% Size of clusters
N_C = accumarray(cluster, 1);

% Check condition to avoid sure singularity of P_gg and M_gg
if n_lev < max(N_C)
    disp('The number of iterations for the estimation of M_gg is too low')
    disp('The estimate of M_gg would be singular.')
    disp(['Proceeding with n_lev = ', num2str(max(N_C) + 10)])
    n_lev = max(N_C) + 10;
end

%Number of elements to be estimated in all the P_gg and M_gg matrices (only
%need to keep lower triangular part to save memory)
nz = sum(N_C.*(N_C + 1)./2);

%Initialize matrices for P_hat_mover and M_hat_mover 
P = spalloc(NT_mover, NT_mover, nz);
M = spalloc(NT_mover, NT_mover, nz);

% Pass variables to workers 
parfor j=1:n_lev
    
    y_lev = randsample([-1;1], NT_mover, true);
    
    % Estimate regression
    [b_lev, ~] = pcg(xx, X'*y_lev, mytol, 300, L);
    
    % Fitted value
    y_hat = X*b_lev;     % fitted values
    
    aux_matrix = sparse(1:NT_mover, cluster, y_hat);
    P = P + tril(aux_matrix * aux_matrix');
    
    aux_matrix = sparse(1:NT_mover, cluster, y_lev - y_hat);    
    M = M + tril(aux_matrix * aux_matrix');   
    
end

% Complete matrices
P = (1./n_lev) .* P;
P_hat = P + P' - diag(diag(P));
clear P

M = (1./n_lev) .* M;
M_hat = M + M' - diag(diag(M));


end