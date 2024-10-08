function [ind_keep, cluster_out] = pruning_cluster_loop(id_orig, firmid_orig, cluster_orig)

%%% To generate final connected set and remove one-timers
firmid_0 = firmid_orig;
firmid = firmid_orig;
id = id_orig;
cluster_0 = cluster_orig;
cluster = cluster_orig;

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
    cluster_0 = cluster_0(move);

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
        firmlst=find(sindex==idx(1)); %firms in connected set
        sel=ismember(firmid,firmlst);
        clear firmlst
        
        % Filter: remaining firmids have a mover by construction
        id = id(sel);
        firmid = firmid(sel);
        cluster = cluster(sel);      
        firmid_0 = firmid_0(sel);
        cluster_0 = cluster_0(sel);
            
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
        cluster_0 = cluster_0(sel);

        % Resetting ids another time
        [~, ~, id] = unique(id);
        [~, ~, firmid] = unique(firmid);    
        [~, ~, cluster] = unique(cluster);         
    end
    
    % Number of unique workers and clusters
    numClusters = max(cluster);
    
    % Check number of elements in bipartite graph 
    size_graph = max(firmid) + max(id);
    
    % Eliminate each cluster at a time and check if subgraph is still
    % connected and has the same elements

    % Preallocate logical array to track bad clusters
    is_bad_cluster = false(numClusters, 1);
    id_par = parallel.pool.Constant(id);
    firmid_par = parallel.pool.Constant(firmid);
    parfor g = 1:numClusters
        % Remove the cluster and recompute the adjacency matrix
        [~, ~, id_aux] = unique(id_par.Value(cluster ~= g));
        [~, ~, firmid_aux] = unique(firmid_par.Value(cluster ~= g));
        A = bipartiteAdjMatrix(id_aux, firmid_aux);
        
        % Compute the size of the largest connected component
        [~, sz_g] = conncomp(graph(A));
        
        % Check if the subgraph is still connected and is the same size
        if ~(size_graph == sz_g(1))
            is_bad_cluster(g) = true; % Mark this cluster as bad
        end
    end

    % Extract indices of bad clusters after the parfor loop
    bad_clusters = find(is_bad_cluster);    
        
    % index of bad clusters
    sel = ~ismember(cluster, bad_clusters);
    n_of_bad_clusters = size(bad_clusters, 1);
    
    % Filter
    id = id(sel);
    firmid = firmid(sel);
    cluster = cluster(sel);
    firmid_0 = firmid_0(sel);
    cluster_0 = cluster_0(sel);

    if isempty(cluster_0)        
        ind_keep = [];
        cluster_out = [];
        return
    end
    
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

%% Observations belonging to bad clusters should get an observation specific cluster
ind_bad_cluster = ~ismember(cluster_orig, cluster_0);

disp(['Number of remaining good clusters: ', num2str(max(cluster))])
disp('Observations belonging to bad clusters get clusters at the observation level')

cluster_out = cluster_orig;
cluster_out(ind_bad_cluster) = (max(cluster_orig) + 1 : max(cluster_orig) + sum(ind_bad_cluster));


