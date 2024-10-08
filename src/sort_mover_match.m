function [ind_mover, cluster, i_sortmover, N_C_long] = sort_mover_match(id,firmid)

    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    % Identify firmid lags
    lagfirmid = [NaN; firmid(1:end-1)];
    lagfirmid(gcs==1) = NaN; % First obs for each worker

    % Identify stayers. Accumulate same firmid cases per id
    stayer = (firmid==lagfirmid);
    stayer(gcs==1) = 1;
    stayer = accumarray(id, stayer);

    % Count observations per id
    T = accumarray(id, 1);
    % Always stayer if same firmid over time for all the id observations
    stayer = T==stayer;
    ind_mover = stayer~=1;
    ind_mover = ind_mover(id);
    
    % Generate clusters at match level
    cluster = findgroups(id, firmid);
    N_C = accumarray(cluster, 1);

    % Bring to NT size (repeat per cluster)
    N_C_long_sort = N_C(cluster);
    N_C_long = N_C(cluster);

    % Assign 0 to stayers
    N_C_long_sort(~ind_mover) = 0;

    % Create cluster information for movers
    [cluster_mover] = findgroups(id(ind_mover), firmid(ind_mover));

    cluster(~ind_mover) = 1; % Assign cluster = 1 for stayers
    cluster(ind_mover) = cluster_mover + 1; % Adjust cluster classification for movers

    % N_C_long_sort to order the variables: Movers in descending order of
    % observations in the match; then stayers
    [~, i_sortmover] = sort(N_C_long_sort, 'descend');

    % Sort ind
    ind_mover   = ind_mover(i_sortmover);
    cluster     = cluster(i_sortmover);

end