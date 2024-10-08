function eps_LC = get_eps_LC_cluster(P, M, r)

    % Do Cholesky factorization
    L = chol(M + P, "lower");
    
    % Solve system to avoid inverse of M
    z = M \ (L *r);
    
    % Get leave-cluster-out errors
    eps_LC = L' * z;

end