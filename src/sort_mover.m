function [ind_mover, i_sortmover] = sort_mover(id,firmid)

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

    % order to put movers first
    [ind_mover, i_sortmover] = sort(ind_mover, 'descend');

end