function [ind_keep] = pruning_match_LdM_old(id, firmid, year)


%%% Keep track of removals
% id_obs will be filtered and renamed
% id_obs_0 will be filtered only
% id_obs_orig will be used to generate the filtering variable of the same size as the original data
NT = size(id,1);
id_obs = [1:NT]';
id_obs_0 = id_obs;

n_of_bad_obs = 1;
ind_drop_LdM = 1;
ind_drop_unique = 1;

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
    id_obs = id_obs(sel);
    year = year(sel);

    % Resetting ids another time
    [~,~,n] = unique(firmid);
    firmid = n;
    [~,~,n] = unique(id);
    id = n;
    
    %%% Remove problematic observations that are not leave-match-out
    %%% estimable
    %find movers
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    move=(firmid~=lagfirmid);
    move(gcs==1)=0;
    move=accumarray(id,move);
    move=(move>0);
    move=move(id);
        
    % Create graph object
    adjMatrix = bipartiteAdjMatrix(id, firmid);
    G = graph(adjMatrix);

    % Dictionary edges and observations
    idxOut = findedge(G, id, firmid + max(id)); % (Nodes corresponding to firms are numbered after nodes corresponding to workers. That's why there is "+ max(id)")

    % Get bridges
    bridges = find_bridges(G);
    
    % Indicator if obs corresponds to bridge AND that obs is a mover.
    % The edges of stayers are bridges; we do not want to remove them
    id_bridge = ismember(idxOut,bridges).* (move == 1);
    n_of_bad_obs = sum(id_bridge);

    % Remove observations that correspond to a bridge or cut-edge
    sel = (id_bridge~=1);

    % filter
    id = id(sel);
    firmid = firmid(sel);
    id_obs = id_obs(sel);
    year = year(sel);

    %%% Remove firm-years with one worker only
    firmid_y = findgroups(firmid, year);
    [GC, GR] = groupcounts([id, firmid_y]);

    dict_id_firmid_weight = [GR{:,1} GR{:,2} GC]; % unique id, unique firmid*year and counts per match*year
    % sort per match*year to repeat later on
    [~, i_sort] = sort(dict_id_firmid_weight(:,2)); % sort per firmid_y for repetition
    dict_id_firmid_weight = dict_id_firmid_weight(i_sort,:);
    
    % Number of observations per firmid*year to decide if LdM moment can be
    % computed and adjust
    N_obs_per_firmid_y = accumarray(dict_id_firmid_weight(:,2),dict_id_firmid_weight(:,3)); % counts of observations per firmid*year
    N_obs_per_firmid_y = N_obs_per_firmid_y(firmid_y); % repeat vector

    ind_drop_LdM = (N_obs_per_firmid_y==1);

    % filter
    id = id(~ind_drop_LdM);
    firmid = firmid(~ind_drop_LdM);
    year = year(~ind_drop_LdM);
    id_obs = id_obs(~ind_drop_LdM);

    %%% Remove ids that appear only once after the pruning. They dont generate connections
    % and have leverage equal to 1
    % Count obs per id
    [~, ~, iunique] = unique(id);
    n = accumarray(iunique(:), 1);

    % Remove observations with only one observation
    ind_drop_unique = n == 1; 
    ind_drop_unique = ind_drop_unique(iunique); % Repeat to filter original vector
    clear iunique n
    
    % Drop one-timers after all the filters
    firmid=firmid(~ind_drop_unique);
    id=id(~ind_drop_unique);
    id_obs = id_obs(~ind_drop_unique);
    year = year(~ind_drop_unique);

    % Rename ids
    [~,~,n]=unique(firmid);
    firmid=n;
    [~,~,n]=unique(id);
    id=n;  

end

% generate indicator to keep of the length of the original data
ind_keep = ismember(id_obs_0,id_obs);
end