function [est_var_r_collapsed, weight_stayers] = collapse_var_stayers_cluster(eps_LC, y, dim_fe_mover, id, cluster, N_C)

    % Rebrand clusters
    cluster = cluster - cluster(1) + 1;
    id = id - dim_fe_mover(1);
    
    NT = sum(N_C>1);

    if NT > 0
        % Create auxiliary matrices
        sp_eps = sparse(1:NT, cluster(1:NT), eps_LC(1:NT));
        sp_y = sparse(1:NT, cluster(1:NT), y(1:NT));
        aux = sp_eps*sp_y';

        % Sum each row. We will later sum each element per cluster and
        % divide by n^2 so we get the variance of the average
        sum_OP = sum(aux, 2);
        
    else
        sum_OP = [];
    end


    sum_OP = full([sum_OP; eps_LC(NT+1:end).*y(NT+1:end)]); %Include values of workers in only one cluster
        
    % Get the number of periods that stayers are in the sample

    NT = length(sum_OP);
    aux = sparse(1:NT,id,1);
    weight_stayers = aux'*ones(NT,1);
    
    est_var_r_collapsed = aux'*sum_OP;
    est_var_r_collapsed = est_var_r_collapsed./(weight_stayers.^2); % Divide by time in sample, squared

    