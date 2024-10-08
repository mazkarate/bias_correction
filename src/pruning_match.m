function ind_keep = pruning_match(id_orig,firmid)


%%% To generate final connected set and remove one-timers
firmid_0 = firmid;
firmid_orig = firmid;
id = id_orig;

n_of_bad_obs = 1;

while n_of_bad_obs >= 1

    % Take largest connected set
    A = build_adj(id,firmid); 
    [sindex, sz] = conncomp(graph(A));
    sz=sz';
    sindex=sindex';
    idx=find(sz==max(sz)); %find largest set
    firmlst=find(sindex==idx); %firms in connected set
    sel=ismember(firmid,firmlst);
    
    % filter: remaining firmids have a mover by construction
    id = id(sel);
    firmid = firmid(sel);
    firmid_0 = firmid_0(sel);

    % Resetting ids another time
    [~, ~, firmid] = unique(firmid);
    [~, ~, id] = unique(id);
        
    % Find movers
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    move=(firmid~=lagfirmid);
    move(gcs==1)=0;
    move=accumarray(id,move);
    move=(move>0);
    move=move(id);
    
    % Work only with movers
    id = id(move);
    firmid = firmid(move);
    firmid_0 = firmid_0(move);

    % Resetting ids of movers another time
    [~, ~, firmid] = unique(firmid);
    [~, ~, id] = unique(id);

    % Create graph object
    adjMatrix = bipartiteAdjMatrix(id, firmid);
    G = graph(adjMatrix);

    % Dictionary edges and observations
    idxOut = findedge(G, id, firmid + max(id)); % (Nodes corresponding to firms are numbered after nodes corresponding to workers. That's why there is "+ max(id)")

    % Get bridges
    bridges = find_bridges(G);
    
    % Indicator if match corresponds to a bridge
    id_bridge = ismember(idxOut, bridges);
    n_of_bad_obs = sum(id_bridge);

    % Remove observations that correspond to a bridge
    sel = (id_bridge~=1);

    % filter
    id = id(sel);
    firmid = firmid(sel);
    firmid_0 = firmid_0(sel);

    % Rename ids
    [~, ~, firmid]=unique(firmid);
    [~, ~, id]=unique(id);    

end

clear id

% Vector to indicate belonging to connected set 
ind_keep = ismember(firmid_orig, firmid_0);

%% Remove one timers after selecting connected set

% Count obs per id. Those workers without a firm in connected set get zeros
n = accumarray(id_orig(ind_keep), 1);

% Find id's that are not one timers after filtering for the connected firms
% Also filter those that are not in the connected set
non_one_timers = find(n > 1);
sel = ismember(id_orig, non_one_timers);

% Get final index: belong to connected set AND is not one timer
ind_keep = ind_keep & sel;


end