function delta = bootstrap_cluster_multi_group(n_boot, B_plus, B_minus, dim_fe_mover, X, xx, L,...
                                        mytol, NT_mover, id, firmid, firmid_stayers, cluster,...
                                        i_group, ch_group, N_u_group)

   N_multi_group = length(N_u_group);
    
    total_length = sum(N_u_group) + N_multi_group;

    % Pre-assign estimate of bias
    delta_v = 0;    

    % Start and end indices
    start_idx = zeros(N_multi_group,1);
    end_idx = zeros(N_multi_group,1);
    for i=1:N_multi_group
        if i == 1
            start_idx(i) = 1;
            end_idx(i) = N_u_group(i) + 1; % To include the whole sample moments
        else
            start_idx(i) = end_idx(i-1) + 1;
            end_idx(i) = start_idx(i) + N_u_group(i);
        end
    end

    NT = size(B_plus,1); %Change length of simulating vectors
    
    % Get number of simulations: 1 per cluster for movers, and one per id for
    % stayers
    n_sim1 = numel(unique(cluster));
    n_sim2 = numel(unique(id(NT_mover + 1 :end)));
    
    parfor s=1:n_boot
        
        % Initialize auxiliary vectors
        aux = zeros(total_length, 3);      

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
        % Clear memory
        b_tilde = []; y_b = [];
        for i = 1:N_multi_group
            pe_tilde = pe_tilde(i_group(:, i));
            fe_tilde = fe_tilde(i_group(:, i));
            aux(start_idx(i):end_idx(i),:) = compute_moments_group(pe_tilde, fe_tilde, ch_group{i});
        end
        
        % B-
        y_b = B_minus.*v_tilde;
        [b_tilde,~]   =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
        [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
        % Clear memory
        b_tilde = []; y_b = []; v_tilde = [];
        for i = 1:N_multi_group
            pe_tilde = pe_tilde(i_group(:, i));
            fe_tilde = fe_tilde(i_group(:, i));
            aux(start_idx(i):end_idx(i),:) = aux(start_idx(i):end_idx(i),:)...
                        - compute_moments_group(pe_tilde, fe_tilde, ch_group{i});                    
        end
    
        % Delta
        delta_v = delta_v + aux;    
       
    end

   % Store bootstrap results
    delta_v = delta_v./n_boot; % take the average

    % Unpack the values of delta into cells
    delta = cell(N_multi_group,1);
    for i = 1:N_multi_group
        delta{i} = delta_v(start_idx(i):end_idx(i),:);
    end

end