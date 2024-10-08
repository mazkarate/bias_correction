function [ind_keep] = pruning_obs_LdM(id, firmid, year_all)


%%% Keep track of removals
% id_obs will be filtered and renamed
% id_obs_0 will be filtered only
% id_obs_orig will be used to generate the filtering variable of the same size as the original data
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


n_of_bad_obs = 1;
ind_drop_LdM = 1;
ind_drop_unique = 1;
j = 0;


while n_of_bad_obs >= 1 || sum(ind_drop_LdM) > 0 || sum(ind_drop_unique) >0 

    %%% Connected set
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
    id_0 = id_0(sel);
	id_obs_mover = id_obs_mover(sel);

    % Filter the sample with movers and stayers
    sel_all = ismember(firmid_all, firmid_0); % indicator of the same length as first variable
    id_all = id_all(sel_all);
    firmid_all = firmid_all(sel_all);
    firmid_y_0 = firmid_y_0(sel_all);
    id_obs_all = id_obs_all(sel_all);

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
	id_obs_mover = id_obs_mover(move);
	if j ==0
		id_obs_stayer = id_obs_all(~move);
	end

    % Resetting ids of movers another time
    [~, ~, firmid] = unique(firmid);
    [~, ~, id] = unique(id);

    %%% Remove problematic observations that are not leave-obs-out
    %%% estimable
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

    % filter sample with movers
    id = id(sel);
    firmid = firmid(sel);
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

     % Rename ids
    [~, ~, firmid]=unique(firmid);
    [~, ~, id]=unique(id); 

	j = j +1;

end

% generate indicator to keep of the length of the original data
ind_keep = ismember(id_obs_orig,id_obs_all);


end