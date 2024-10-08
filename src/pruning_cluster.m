function ind_keep = pruning_cluster(id_orig, firmid, cluster)

%%% To generate final connected set and remove one-timers
firmid_0 = firmid;
firmid_orig = firmid;
id = id_orig;

n_of_bad_clusters = 1;

while n_of_bad_clusters >= 1
    %%% Remove workers and firms that are contained in just one cluster
    n_of_bad_w_f = 1;

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
    cluster = cluster(move);
    firmid_0 = firmid_0(move);

    % Resetting ids of movers 
    [~, ~, firmid] = unique(firmid);
    [~, ~, id] = unique(id);
    [~, ~, cluster] = unique(cluster);
   
    while n_of_bad_w_f >= 1
       
        % Take largest connected set
        A = build_adj(id, firmid); 
        [sindex, sz] = conncomp(graph(A));
        sz=sz';
        sindex=sindex';
        idx=find(sz==max(sz)); %find largest set
        firmlst=find(sindex==idx); %firms in connected set
        sel=ismember(firmid,firmlst);
        clear firmlst
        
        % Filter: remaining firmids have a mover by construction
        id = id(sel);
        firmid = firmid(sel);
        cluster = cluster(sel);      
        firmid_0 = firmid_0(sel);
            
        % Resetting ids
        [~, ~, id] = unique(id);
        [~, ~, firmid] = unique(firmid);    
        [~, ~, cluster] = unique(cluster);        

        % Create graph object
        adjMatrix = tripartiteAdjMatrix(cluster, id, firmid);
        G = graph(adjMatrix);
        
        % Compute the degree of each node
        nodeDegrees = degree(G);
        
        % Number of unique workers and clusters
        numClusters = max(cluster);
        numWorkers = max(id);
        
        % Select degress for workers and firms nodes
        nodeDegrees_id = nodeDegrees(numClusters + 1 : numClusters + numWorkers);
        bad_id = find(nodeDegrees_id == 1);
        nodeDegrees_firmid = nodeDegrees(numClusters + numWorkers + 1 : end);
        bad_firmid = find(nodeDegrees_firmid == 1);
        
        % Select bad id and bad firmid 
        sel_bad_id = ismember(id, bad_id); 
        sel_bad_firmid = ismember(firmid, bad_firmid); 
        n_of_bad_w_f = sum(sel_bad_id) + sum(sel_bad_firmid);

        % Select observations that do not belong to bad workers or bad
        % firms
        sel = ~(sel_bad_id | sel_bad_firmid);
        
        % Filter
        id = id(sel);
        firmid = firmid(sel);
        cluster = cluster(sel);        
        firmid_0 = firmid_0(sel);        

        % Resetting ids another time
        [~, ~, id] = unique(id);
        [~, ~, firmid] = unique(firmid);    
        [~, ~, cluster] = unique(cluster);         
    end
    
    % Number of unique workers and clusters
    numClusters = max(cluster);    

    % Create graph object after deleting bad workers 
    adjMatrix = tripartiteAdjMatrix(cluster, id, firmid);    

    % Get cluster nodes that are articulation points
    [~, artic_points] = biconncomp(graph(adjMatrix)); % Matlab built in function
    bad_cluster = artic_points(artic_points <= numClusters)';
    sel = ~ismember(cluster, bad_cluster);
    n_of_bad_clusters = size(bad_cluster, 1);
    
    % Filter
    id = id(sel);
    firmid = firmid(sel);
    cluster = cluster(sel);
    firmid_0 = firmid_0(sel);
    
    % Resetting ids one last time
    [~, ~, id] = unique(id);
    [~, ~, firmid] = unique(firmid);    
    [~, ~, cluster] = unique(cluster); 
    

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