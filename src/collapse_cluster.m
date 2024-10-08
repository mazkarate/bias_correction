function [y_col, id_col, firmid_col, weights, ind_mover_col] = collapse_cluster(y, id, firmid, cluster, dimensions_mover)

    % Indicator of stayer (last ids) 
    ind_stayer = id>dimensions_mover(1);

    % Create cluster indicator for only movers
    cluster_norm = cluster(~ind_stayer)-1; % Clusters of stayers is equal to 1    

    %Average outcome variable for movers
    y_col = accumarray(cluster_norm, y(~ind_stayer), [], @(x) mean(x));
    
    %Get length of vector of movers after collapsing at the cluster level
    N_mover_col = size(y_col,1);

    %Get id and firmid values for each cluster
    id_col = accumarray(cluster_norm, id(~ind_stayer), [], @(x) x(1));
    firmid_col = accumarray(cluster_norm, firmid(~ind_stayer), [], @(x) x(1));
    
    %Get number of ocurrences per cluster
    weights = accumarray(cluster_norm, 1);

    %Add values for stayers
    y_col = [y_col; y(ind_stayer)];
    id_col = [id_col; id(ind_stayer)];
    firmid_col = [firmid_col; firmid(ind_stayer)];
    weights = [weights; ones(sum(ind_stayer), 1)]; % Add weight of one for stayers
    ind_mover_col = logical([ones(N_mover_col,1); zeros(sum(ind_stayer),1)]); % Change ind_movers after collapse per cluster