function ind_keep = pruning_obs_LdM(id_orig, firmid, year_all)


%%% To generate final connected set and remove one-timers
firmid_0 = firmid;
firmid_orig = firmid;
id = id_orig;
id_0 = id_orig;
firmid_all = firmid;
id_all = id;
firmid_y_orig = findgroups(firmid, year_all);
firmid_y_0 = firmid_y_orig;

n_of_bad_obs = 1;

while n_of_bad_obs >= 1

    %%% Remove firm-years with one worker only
    [GC, GR] = groupcounts([id_all, firmid_y_0]);

    dict_id_firmid_weight = [GR{:,1} GR{:,2} GC]; % unique id, unique firmid*year and counts per match*year
    % % sort per match*year to repeat later on
    % [~, i_sort] = sort(dict_id_firmid_weight(:,2)); % sort per firmid_y for repetition
    % dict_id_firmid_weight = dict_id_firmid_weight(i_sort,:);
    
    % Number of observations per firmid*year to decide if LdM moment can be
    % computed and adjust
    N_obs_per_firmid_y = accumarray(dict_id_firmid_weight(:,2),dict_id_firmid_weight(:,3)); % counts of observations per firmid*year
    N_obs_per_firmid_y = N_obs_per_firmid_y(firmid_y_0); % repeat vector

    ind_drop_LdM = (N_obs_per_firmid_y==1);

    % filter the sample with movers and stayers
    id_all = id_all(~ind_drop_LdM);
    firmid_all = firmid_all(~ind_drop_LdM);
    year_all = year_all(~ind_drop_LdM);
    firmid_y_0 = firmid_y_0(~ind_drop_LdM);

    % filter the sample with movers only
    sel_f = ismember(firmid_0, firmid_all); % indicator of the same length as firmid
    sel_w = ismember(id_0, id_all); % indicator of the same length as id
    sel = sel_f & sel_w;

    id = id(sel);
    firmid = firmid(sel);
    firmid_0 = firmid_0(sel);
    id_0 = id_0(sel);

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
    % Sample with movers and stayers
    firmid_all = firmid_all(~ind_drop_unique);
    id_all = id_all(~ind_drop_unique);
    year_all = year_all(~ind_drop_unique);
    firmid_y_0 = firmid_y_0(~ind_drop_unique);

    % filter the sample with movers only
    sel_f = ismember(firmid_0, firmid_all); % indicator of the same length as firmid
    sel_w = ismember(id_0, id_all); % indicator of the same length as id
    sel = sel_f & sel_w;

    id = id(sel);
    firmid = firmid(sel);
    firmid_0 = firmid_0(sel);
    id_0 = id_0(sel);

    %%% Take largest connected set
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
    id_0 = id_0(sel);

     % Filter the sample with movers and stayers
    sel_all = ismember(firmid_all, firmid_0); % indicator of the same length as first variable
    id_all = id_all(sel_all);
    firmid_all = firmid_all(sel_all);
    year_all = year_all(sel_all);
    firmid_y_0 = firmid_y_0(sel_all);

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
    id_0 = id_0(move);


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
    
    % Indicator if obs corresponds to bridge
    id_bridge = ismember(idxOut,bridges);
    
    % Get weights
    weights = G.Edges.Weight(idxOut);

    % Indicator leverage equal to one
    ind_lev_unit = id_bridge .* (weights == 1);
    n_of_bad_obs = sum(ind_lev_unit);

    % Remove observations with leverage equal to one
    sel = (ind_lev_unit~=1);

    % filter the sample with movers only
    id = id(sel);
    firmid = firmid(sel);
    firmid_0 = firmid_0(sel);
    id_0 = id_0(sel);


    % filter the sample with movers and stayers
    sel_f = ismember(firmid_all, firmid_0); % indicator of the same length as first variable
    sel_w = ismember(id_all, id_0); 
    sel_all = sel_f & sel_w;

    id_all = id_all(sel_all);
    firmid_all = firmid_all(sel_all);
    year_all = year_all(sel_all);
    firmid_y_0 = firmid_y_0(sel_all);

    disp('Size id_all: ',num2str(length(id_all)))

    % Rename ids
    [~, ~, firmid]=unique(firmid);
    [~, ~, id]=unique(id); 

end


% Vector to indicate belonging to connected set 
ind_keep = ismember(firmid_orig, firmid_0);

%% Remove one timers and make LdM restriction after selecting connected set
%%% Loop to filter workers and observations. Make sure that we don't have one timers nor unique firm*year
%%% observations
id_con = id_orig(ind_keep);
firmid_y_con = firmid_y_orig(ind_keep);

disp('Size id_con: ',num2str(length(id_con)))

% firmid_con = firmid_orig(ind_keep);
% 
% % number of firmids and workers
% n_firmid = length(unique(firmid_con));
% n_worker = length(unique(id_con));
% disp(['Number of unique firmids in connected: ', num2str(n_firmid)])
% disp(['Number of unique workers in connected: ', num2str(n_worker)])

diff_size = 1;

while diff_size > 0
    % Count obs per id. Those workers without a firm in connected set get zeros
    n = accumarray(id_con, 1);
    
    % Find id's that are not one timers after filtering for the connected firms
    % Also filter those that are not in the connected set
    non_one_timers = find(n > 1);
    sel = ismember(id_con, non_one_timers);

    size_one_timers = sum(sel);


    % update vectors
    firmid_y_con = firmid_y_con(sel);
    id_con = id_con(sel);
    
    
    %%% Remove firm-years with one worker only
    [GC, GR] = groupcounts([id_con, firmid_y_con]);
    
    dict_id_firmid_weight = [GR{:,1} GR{:,2} GC]; % unique id, unique firmid*year and counts per match*year
    % sort per match*year to repeat later on
    [~, i_sort] = sort(dict_id_firmid_weight(:,2)); % sort per firmid_y for repetition
    dict_id_firmid_weight = dict_id_firmid_weight(i_sort,:);
    
    % Number of observations per firmid*year to decide if LdM moment can be
    % computed and adjust
    N_obs_per_firmid_y = accumarray(dict_id_firmid_weight(:,2),dict_id_firmid_weight(:,3)); % counts of observations per firmid*year
    N_obs_per_firmid_y = N_obs_per_firmid_y(firmid_y_con); % repeat vector
    
    ind_keep_LdM = (N_obs_per_firmid_y>=2);
    
    % filter of firm*year observations
    firmid_y_filt = firmid_y_con(ind_keep_LdM);
    
    % create filter of original sample size
    ind_keep_LdM = ismember(firmid_y_con, firmid_y_filt);

    size_LdM = sum(ind_keep_LdM);

    % update vectors
    firmid_y_con = firmid_y_con(ind_keep_LdM);
    id_con = id_con(ind_keep_LdM);

    diff_size = abs(size_one_timers - size_LdM);

end

%%% Worker and observation filter
sel = ismember(id_orig, id_con);
ind_keep_LdM = ismember(firmid_y_orig, firmid_y_con);

%%% Get final index: belong to connected set AND is not one timer AND
%%% survive LdM filters
ind_keep = ind_keep & sel & ind_keep_LdM;

disp('Size final sample: ',num2str(length(id_orig(ind_keep))))


end