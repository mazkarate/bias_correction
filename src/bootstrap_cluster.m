function [delta, corrected] = bootstrap_cluster(n_boot, B_plus, B_minus, dim_fe_mover, X, xx, L,...
                                        plugin, mytol, NT_mover, id, firmid, firmid_stayers, cluster)

    % Pre-assign estimate of bias
    delta = 0;
    
    % Get number of simulations: 1 per cluster for movers, and one per id for
    % stayers
    n_sim1 = numel(unique(cluster));
    n_sim2 = numel(unique(id(NT_mover + 1 :end)));
    
    parfor s=1:n_boot
        
        % Determine the dependent variable: first cluster 
        v_tilde = randsample([-1;1],n_sim1,true);
        v_tilde2 = randsample([-1;1],n_sim2,true);

        % Repeat values for movers for the different clusters and stack the
        % rest
        v_tilde = [v_tilde(cluster); v_tilde2];
    
        % Solve normal equations with preconditioning matrix L
        % B+
        y_b = B_plus.*v_tilde;
        [b_tilde,~]    =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
        [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
        aux = compute_moments(pe_tilde, fe_tilde);
        
        % B-
        y_b = B_minus.*v_tilde;
        [b_tilde,~]   =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
        [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
        aux = aux - compute_moments(pe_tilde, fe_tilde);
        % Delta                        
        delta = delta + aux;
       
    end

    %Store bootstrap results
    delta = delta./n_boot; % take the average
	corrected = plugin - delta;

end