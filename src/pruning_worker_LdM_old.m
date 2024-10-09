function [ind_keep] = pruning_worker_LdM_old(id, firmid, year)
%%%% Adapted from the replication package of the Kline, P., Saggio, R., & Sølvsten, M. (2020). Leave-out estimation of variance components. Econometrica, 88(5), 1859-1898.
%%% Source: https://www.dropbox.com/s/iaj3aap3ibfhup8/Replication%20ECTA.zip?dl=1 function pruning_unbal_v3.m


%  This function finds the leave one out largest connected set. 
%  That is, a connected set where if we were to remove the entire working
%  history of a given individual, the resulting graph still remains
%  connected. See appendix B of KSS.
%     
%-------- INPUTS: 
% id: worker identifiers. Dimensions: N* x 1
% firmid: firm identifiers. Dimensions: N* x 1

%  -------- OUTPUT:    
% ind_keep: indicator to keep observation that will filter the original data. will need
% to update the identifier codes outside

% %%% Keep track of removals
% % id will be filtered and renamed
% % id_0 will be filtered only
% % id_orig will be used to generate the filtering variable of the same size as the original data
% id_0    = id;
% id_orig = id_0; % to come back to the original data size

%%% Keep track of removals
% id_obs will be filtered and renamed
% id_obs_0 will be filtered only
% id_obs_orig will be used to generate the filtering variable of the same size as the original data
NT = size(id,1);
id_obs = [1:NT]';
id_obs_0 = id_obs;

n_of_bad_workers = 1;
ind_drop_LdM = 1;
ind_drop_unique = 1;

while n_of_bad_workers >= 1 || sum(ind_drop_LdM) > 0 || sum(ind_drop_unique) >0

    %%% Remove workers that are problematic
    %find movers
    J=max(firmid);
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    move=(firmid~=lagfirmid);
    move(gcs==1)=0;
    move=accumarray(id,move);
    move=(move>0);
    move=move(id);
    
    %id and firmids associated with movers only
    id_movers=id(move);
    firmid_mover=firmid(move);

    %need to normalize id_movers and keep a dictionary. 
    id_movers_orig=id_movers;
    [~,~,id_movers]=unique(id_movers);
    n_movers=max(id_movers);
    
    %unique pairs    
    list=[id_movers firmid_mover];
    list=unique(list,'rows');
    
    %Construct Adj. Matrix of Bipartite Graph
    B=sparse(list(:,1),list(:,2),1,n_movers,J);
    G = [ sparse( n_movers, n_movers ), B; B.', sparse(J, J)];
    
    %Inspect the resulting Graph: find the workers that constitute an
    %articulation point
    [~,artic_points] = biconncomp(graph(G)); % Matlab built in function
    bad_workers=artic_points(artic_points<=n_movers);
    
    %Now get the right index for these workers
    sel=ismember(id_movers,bad_workers);
    bad_workers=id_movers_orig(sel);
    bad_workers=unique(bad_workers);
    n_of_bad_workers=size(bad_workers,1);

    %Drop these workers
    sel=~ismember(id,bad_workers);
    id=id(sel);
    firmid=firmid(sel);
    % id_0 = id_0(sel);
    id_obs = id_obs(sel);
    year = year(sel);
    
	% Rename ids
	[~,~,n]=unique(firmid);
	firmid=n;
	[~,~,n]=unique(id);
	id=n;
 
	% Adjacency matrix
    A=build_adj(id,firmid);
    
    [sindex,sz] = conncomp(graph(A)); % Matlab built in function. Get connected sets
    idx=find(sz==max(sz)); %find largest set
    firmlst=find(sindex==idx); %firms in connected set
    sel=ismember(firmid,firmlst);
    
	% Select the connected set
    firmid=firmid(sel);
    id=id(sel);
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
    [~,~,n] = unique(firmid);
    firmid = n;
    [~,~,n] = unique(id);
    id=n;  
 
end

% generate indicator to keep of the length of the original data
ind_keep = ismember(id_obs_0,id_obs); % search A in B

end

    