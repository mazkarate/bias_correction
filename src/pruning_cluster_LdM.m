function ind_keep = pruning_cluster_LdM(id_orig, firmid, cluster, year_all)

%%% To generate final connected set and remove one-timers
id = id_orig;
NT = size(id,1);
id_obs_all = [1:NT]';
id_obs_orig = id_obs_all;
id_obs_mover = id_obs_orig;
firmid_0 = firmid;

id_0 = id;
firmid_all = firmid;
id_all = id;
firmid_y_orig = findgroups(firmid, year_all);
firmid_y_0 = firmid_y_orig;


n_of_bad_clusters = 1;
ind_drop_LdM = 1;
ind_drop_unique = 1;
j = 0;

while n_of_bad_clusters >= 1 || sum(ind_drop_LdM) > 0 || sum(ind_drop_unique) >0 
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
    id_0 = id_0(move);
	id_obs_mover = id_obs_mover(move);
	if j ==0
		id_obs_stayer = id_obs_all(~move);
    end
    % mover = move(move);


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
        id_0 = id_0(sel);
	    id_obs_mover = id_obs_mover(sel);
    
        % Filter the sample with movers and stayers
        sel_all = ismember(firmid_all, firmid_0); % indicator of the same length as first variable
        id_all = id_all(sel_all);
        firmid_all = firmid_all(sel_all);
        firmid_y_0 = firmid_y_0(sel_all);
        id_obs_all = id_obs_all(sel_all);
            
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
        id_0 = id_0(sel);
	    id_obs_mover = id_obs_mover(sel);
    
        %%% filter the sample with movers and stayers
	    % create filtered id_obs
	    sel_mover = ismember(id_obs_all, id_obs_mover); % movers
	    sel_stayer = ismember(id_obs_all, id_obs_stayer); % stayer
	    
	    sel_all = sel_mover + sel_stayer;
        sel_all = logical(sel_all);
	    
        id_all = id_all(sel_all);
        firmid_all = firmid_all(sel_all);
        firmid_y_0 = firmid_y_0(sel_all);
        id_obs_all = id_obs_all(sel_all);

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
    
    % Filter sample with movers
    id = id(sel);
    firmid = firmid(sel);
    cluster = cluster(sel);
    firmid_0 = firmid_0(sel);
    id_0 = id_0(sel);
	id_obs_mover = id_obs_mover(sel);

    %%% filter the sample with movers and stayers
	% create filtered id_obs
	sel_mover = ismember(id_obs_all, id_obs_mover); % movers
	sel_stayer = ismember(id_obs_all, id_obs_stayer); % stayer
	
	sel_all = sel_mover + sel_stayer;
    sel_all = logical(sel_all);
	
    id_all = id_all(sel_all);
    firmid_all = firmid_all(sel_all);
    firmid_y_0 = firmid_y_0(sel_all);
    id_obs_all = id_obs_all(sel_all);

    %%% Remove firm-years with one worker only
    [GC, GR] = groupcounts([id_all, firmid_y_0]);

    dict_id_firmid_weight = [GR{:,1} GR{:,2} GC]; % unique id, unique firmid*year and counts per match*year
    % sort per match*year to repeat later on
    [~, i_sort] = sort(dict_id_firmid_weight(:,2)); % sort per firmid_y for repetition
    dict_id_firmid_weight = dict_id_firmid_weight(i_sort,:);
    
    % Number of observations per firmid*year to decide if LdM moment can be
    % computed and adjust
    N_obs_per_firmid_y = accumarray(dict_id_firmid_weight(:,2),dict_id_firmid_weight(:,3)); % counts of observations per firmid*year
    N_obs_per_firmid_y = N_obs_per_firmid_y(firmid_y_0); % repeat vector

    ind_drop_LdM = (N_obs_per_firmid_y==1);

    % filter the problematic for LdM. Movers and stayers
    id_all  = id_all(~ind_drop_LdM);
    firmid_all = firmid_all(~ind_drop_LdM);
    firmid_y_0 = firmid_y_0(~ind_drop_LdM);
    id_obs_all = id_obs_all(~ind_drop_LdM);

	%%% filter the sample with movers
	% create filtered id_obs
	sel = ismember(id_obs_mover, id_obs_all); % movers

    % filter sample with movers
    id = id(sel);
    firmid = firmid(sel);
    firmid_0 = firmid_0(sel);
    id_0 = id_0(sel);
	id_obs_mover = id_obs_mover(sel);

	%%% Filter observations of stayers
	sel_stayer = ismember(id_obs_stayer, id_obs_all); % stayer
	id_obs_stayer = id_obs_stayer(sel_stayer);

    %%% Remove ids that appear only once after the pruning. They dont generate connections
    % and have leverage equal to 1
    % Count obs per id
    [~, ~, iunique] = unique(id_all);
    n = accumarray(iunique(:), 1);

    % Remove observations with only one observation
    ind_drop_unique = n == 1; 
    ind_drop_unique = ind_drop_unique(iunique); % Repeat to filter original vector
    clear iunique n
    
    % Drop one-timers after all the filters
    firmid_all = firmid_all(~ind_drop_unique);
    id_all = id_all(~ind_drop_unique);
    firmid_y_0 = firmid_y_0(~ind_drop_unique);
    id_obs_all = id_obs_all(~ind_drop_unique);

    % filter the sample with movers only
    sel = ismember(id_obs_mover, id_obs_all); % indicator of the same length as firmid

    id = id(sel);
    firmid = firmid(sel);
    firmid_0 = firmid_0(sel);
    id_0 = id_0(sel);
	id_obs_mover = id_obs_mover(sel);
	
	%%% Filter observations of stayers
	sel_stayer = ismember(id_obs_stayer, id_obs_all); % stayer
	id_obs_stayer = id_obs_stayer(sel_stayer);
   
    % Resetting ids one last time
    [~, ~, id] = unique(id);
    [~, ~, firmid] = unique(firmid);    
    [~, ~, cluster] = unique(cluster); 

	j = j +1;
   
end

% generate indicator to keep of the length of the original data
ind_keep = ismember(id_obs_orig,id_obs_all);

end