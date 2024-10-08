function A=build_adj(id,firmid)
%%%% Adapted from code from Kline, P., Saggio, R., & Sølvsten, M. (2020). Leave‐out estimation of variance components. Econometrica, 88(5), 1859-1898.
%%%% Source: https://github.com/rsaggio87/LeaveOutTwoWay/blob/master/codes/build_adj.m


%%% Builds adjacency matrix taking the movers only

    % Index to find first observation per id
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    % Identify firmid lags
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    % Identify stayers. Accumulate same firmid cases per id
    stayer=(firmid==lagfirmid);
    stayer(gcs==1)=1;
    stayer=accumarray(id,stayer);
    % Count observations per id
    T=accumarray(id,1);
    % Always stayer if same firmid over time for all the id observations
    stayer=T==stayer;
    movers=stayer~=1;
    movers=movers(id);
    % Get ids of movers
    id_movers=id(movers);
    id_movers_orig=id_movers;
    %[ids,m,id_movers]=unique(id_movers);
    % Get firmids of movers, firms potentially in the main connected set
    firmid_movers=firmid(movers);
    % Index to find first observation per id of the movers
    gcs_mover=gcs(movers);
    lagfirmid=[NaN; firmid_movers(1:end-1)];
    lagfirmid(gcs_mover==1)=NaN; %%first obs for each worker
    % Main list of movers
    list=[lagfirmid firmid_movers id_movers_orig]; %%all the moves observed from one period to the next this time only for movers
    sel=~isnan(list(:,1));
    % Main list of movers taking out first observations per id
    list=list(sel,:);
    % Main list of moves from the movers
    sel=list(:,1)~=list(:,2);
    list=list(sel,:);
    % Get unique identifiers to reduce the dimension of vectors forming the sparse adjacency matrix
    [u_list,i2unique,i2orig]=unique(list(:,1:2),'rows');
    weight=histc(i2orig,1:numel(i2unique)); %frequencies
    J=max(firmid);
    % Adjacency matrix
    A = sparse([u_list(:,1);u_list(:,2)],[u_list(:,2);u_list(:,1)],[weight;weight],J,J); %adjacency matrix
end

