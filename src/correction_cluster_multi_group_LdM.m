function [plugin, delta, corrected, decomp_pi, decomp_b, dim_fe, NT] = correction_cluster_multi_group_LdM(y, id, firmid, cluster, year, group, other_fe, mat_controls, n_lev,...    
                                                                                                n_boot, ind_light, mytol, ind_export, filename, v_filename_group)
    %%%% CORRECTION FOR CLUSTERS GIVEN BY THE USER %%%%%
    % Description: It only does the bias corrections on the entire sample.
    % Suited to do leave-cluster-out variance estimate.
    % For corrections for different groups of the sample use "correction_cluster_group.m"
    
    % ----------------------------------------------------------------------
    %   START
    % ----------------------------------------------------------------------
    
    
    disp('Setup: we use the leave-cluster-out variance estimate, so covariance matrix is block diagonal')

    %% Initial variables
    disp(['Size initial data: ',num2str(length(id))])
    
    % Number of FE
    n_FE = 2 + size(other_fe, 2);

    % Store group
    group_old = group;

    % Lagfirmid
    gcs = [NaN; double(id(1:end-1))];
    gcs = id~=gcs;
    lagfirmid = [NaN; double(firmid(1:end-1))];
    lagfirmid(gcs==1) = NaN; %%first obs for each worker
    clear gcs

    % Rename FE identifiers
    disp('Relabeling ids...')
    NT = length(y);
    sel = ~isnan(lagfirmid);

    % Relabel the firms
    [~, ~, n] = unique([firmid; lagfirmid(sel)]);

    firmid = n(1:NT);
    lagfirmid(sel) = n(NT+1:end);
    
    % Relabel the workers
    [~, ~, id] = unique(id);      

	%% Compute the connected set and relabel
    disp('---------- Connected set ---------- ')

    % Find the connected set. Then define dimensions and relevant matrices.
    sel = connected(firmid, lagfirmid, sel);
    clear lagfirmid

    % Filter after finding connected set
    y=y(sel); firmid=firmid(sel); id=id(sel); year=year(sel);
    cluster = cluster(sel);

    if ~isempty(mat_controls)
        mat_controls = mat_controls(sel,:); 
        if ind_light ==1
            % Clear memory from base and caller workspace
            evalin('base', 'mat_controls = [];'); evalin('base', 'mat_controls = [];')
        end
    end    
    if n_FE > 2
        other_fe = other_fe(sel,:);
        if ind_light ==1
            % Clear memory from base and caller workspace
            evalin('base', 'other_fe = [];'); evalin('base', 'other_fe = [];')
        end
    end
    group = group(sel,:); 
    group_old = group_old(sel,:);

    % Relabel ids after keeping the largest connected set
    disp('Relabeling ids again...')           
    [~, ~, firmid] = unique(firmid);     
    [~, ~, id] = unique(id);
    [~, ~, cluster] = unique(cluster);    

    N_multi_group = size(group,2);

    if ind_light ==1
        % Clear memory from base and caller workspace
        evalin('base', 'id = [];'); evalin('base', 'id = [];')
        evalin('base', 'firmid = [];'); evalin('base', 'firmid = [];')
        evalin('base', 'cluster = [];'); evalin('base', 'cluster = [];')
        evalin('base', 'group = [];'); evalin('base', 'group = [];')
        evalin('base', 'y = [];'); evalin('base', 'y = [];')
    end

    disp(['Size connected set data: ',num2str(size(y,1))])
        
    %% Pruning of data and removal of one-time observations
    
    disp('Type of data selection procedure: leave-cluster-out')
    
    fprintf('\n')
    disp('---------- Pruning ---------- ')
    tic
    
    % Check what type of pruning to do: using tripartite graph or looping
    % over the bipartite graph
    firm_worker_cluster = check_type_cluster(id, firmid, cluster);

    % Do pruning at the cluster level: selects surviving firms and clusters
    if firm_worker_cluster
        sel = pruning_cluster_LdM(id, firmid, cluster, year);
    else
        % New cluster vector gets cluster at observation level for those
        % belonging to bad clusters
        [sel, cluster] = pruning_cluster_loop_LdM(id, firmid, cluster, year);   
        if isempty(sel)
            fprintf('\n')
            error('ERROR: Sample is not leave-cluster-out estimable. Consider different level of clustering.' );            
        end
    end
    
    % Filter variables and rename identifiers again
    y = y(sel); id = id(sel); firmid = firmid(sel); year=year(sel);
    cluster = cluster(sel);
    if ~isempty(mat_controls)
        mat_controls = mat_controls(sel, :);   
    end
    if n_FE > 2
        other_fe = other_fe(sel, :);
    end
    group = group(sel,:);
    group_old = group_old(sel,:);
    
    % Rename identifiers
    [~, ~, firmid] = unique(firmid, 'stable');
    [~, ~, id] = unique(id, 'stable');
    [~, ~, cluster] = unique(cluster,  'stable');
    if n_FE > 2
        for j = 1 : n_FE - 2
            [~, ~ ,n] = unique( other_fe(:, j), 'stable' );
            other_fe(:, j) = n;
        end
    end

    toc        
    disp(['Size of pruned data: ', num2str(size(y,1))])
    
    clear sel n    
    
    %% Normalized dimensions, design matrices. Movers and stayers
    NT = length(id);    

    % Raw and normalized dimensions of FE
    dim_fe = [max(id) max(firmid)-1];    
    
    % Initial design matrix with full final sample (movers and stayers)
    % Design matrix first two fixed effects
    D   = sparse(1:NT, id', 1);
    F   = sparse(1:NT, firmid', 1);
	
    % Need to normalize one firm effect effect to avoid collinearity
    % Take out one firm effect. Last firm
    F(:, end) = [];    

    %% Build matrix of dummies for the other FE
    if n_FE > 2
        disp('Building matrix of dummies for the additional fixed effects...')
        C = []; % initialize matrix of dummies for all fixed effects
        for j = 1 : n_FE - 2
            tmp = sparse(1:NT, double(other_fe(:,j))', 1);
            tmp = tmp(:,1 : end-1); % Take out one fixed effect to avoid collinearity
            C = [C, tmp];
        end
        mat_controls = [mat_controls, C];
        clear C other_fe tmp
        disp('Done!')
    end

	%% Residualize extra Fixed Effects and Controls 
    
    if exist('mat_controls', 'var')  && isempty(mat_controls) == 0
       disp('Residualizing controls and/or extra fixed effects')
       
       X  = [D, F, mat_controls]; % Movers and stayers
       xx = X'*X;
       xy = X'*y;
        
       % Degrees of freedom
       dof = NT- size(X,2);

       % Preconditioner for general matrices
       L =lchol_iter(xx);
       if size(L,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
          % Estimation with all covariates
          b = pcg(xx,xy,mytol,300,L,L');
       else 
          b = pcg(xx,xy,mytol,300);
       end

       % Redefine left-hand variable. Variance decompositions will be done 
       % on this variable
       n_2FE = sum(dim_fe); % Dimension of two leading fixed effects
       y = y - mat_controls*b(n_2FE+1:end);
       clear mat_controls

       disp('Done!')
    else          
       disp('No controls nor extra fixed effects to residualize...')
       disp('Proceeding with two leading fixed effects')

       % Degrees of freedom
       dof = NT- sum(dim_fe);       
    end

    %% Define stayers and movers to work with them separately
    % Sort movers at the beginning of the sample
    [ind_mover, i_sortmover] = sort_mover(id, firmid); % ind_mover is already sorted (movers come first)
            
    % Order the rest of the variables
    y = y(i_sortmover); id = id(i_sortmover);
    firmid = firmid(i_sortmover); cluster = cluster(i_sortmover); year = year(i_sortmover);
    group = group(i_sortmover,:); group_old = group_old(i_sortmover,:);
    clear i_sortmover  

    % Rename identifiers again
    [~, ~, firmid] = unique(firmid, 'stable');
    [~, ~, id] = unique(id, 'stable');
    [~, ~, cluster] = unique(cluster, 'stable');

    
    % Get how many movers are
    NT_mover    = sum(ind_mover);
    dim_fe_mover    = [max(id(ind_mover)) max(firmid(ind_mover))-1];
    disp(['Share of mover ids: ', num2str(dim_fe_mover(1)/dim_fe(1))])

    %%% Get unique ids and firmid
    N_id = numel(unique(id));
    N_firmid = numel(unique(firmid));

    %% Define workers that belong to one cluster. Use leave-observation-out variance estimate with them
    % Sort workers with multiple clusters at the beginning of the sample
    % for movers and then for stayers
    [ind_mover_mult_cluster, i_sort_mover_mult_cluster] = sort_mover(id(ind_mover), cluster(ind_mover));

    [ind_stayer_mult_cluster, i_sort_stayer_mult_cluster] = sort_mover(id(~ind_mover), cluster(~ind_mover));

    i_sort_mult_cluster = [i_sort_mover_mult_cluster; i_sort_stayer_mult_cluster + length(i_sort_mover_mult_cluster)];    
    clear i_sort_mover_mult_cluster i_sort_stayer_mult_cluster

    % Order the rest of the variables
    y = y(i_sort_mult_cluster); id = id(i_sort_mult_cluster); year = year(i_sort_mult_cluster);
    firmid = firmid(i_sort_mult_cluster); cluster = cluster(i_sort_mult_cluster);
    clear i_sort_mult_cluster  

    % Rename identifiers another time
    [~, ~, firmid] = unique(firmid, 'stable');
    [~, ~, id] = unique(id, 'stable');
    [u_cluster, ~, cluster] = unique(cluster, 'stable');

    N_clusters = numel(u_cluster);
    clear u_cluster

    ind_mult_cluster = [ind_mover_mult_cluster; ind_stayer_mult_cluster];

    % Get how many movers and stayers have multiple clusters
    NT_mover_mult_cluster    = sum(ind_mover_mult_cluster);    
    dim_fe_mover_mult_cluster    = [length(unique((id(ind_mult_cluster & ind_mover)))), length(unique((firmid(ind_mult_cluster & ind_mover))))-1];
    if NT_mover_mult_cluster > 0
        disp(['Share of mover ids with multiple clusters out of whole sample: ', num2str(dim_fe_mover_mult_cluster(1)/dim_fe(1))])
    else
        disp('Share of mover ids with multiple clusters out of whole sample: 0')
    end


    NT_stayer_mult_cluster    = sum(ind_stayer_mult_cluster);
    dim_fe_stayer_mult_cluster    = [length(unique((id(ind_mult_cluster & ~ind_mover)))), length(unique((firmid(ind_mult_cluster & ~ind_mover))))-1];
    if NT_stayer_mult_cluster > 0
        disp(['Share of stayer ids with multiple clusters out of whole sample: ', num2str(dim_fe_stayer_mult_cluster(1)/dim_fe(1))])
    else
        disp('Share of stayer ids with multiple clusters out of whole sample: 0')
    end


    %% Sort movers and stayers by cluster size.
    % Accomodate data to be ordered also by clusters to get block diagonal
    % matrices

    % Observations of Workers that belong to one cluster become an
    % individual cluster (cluster per observation) this means we are doing
    % leave-observation-out for those
    cluster(NT_mover_mult_cluster + 1 : NT_mover) = (N_clusters + 1 : N_clusters + sum(~ind_mover_mult_cluster))';
    N_clusters = N_clusters + sum(~ind_mover_mult_cluster);
    cluster(NT_mover + NT_stayer_mult_cluster + 1 : NT) = (N_clusters + 1 : N_clusters + sum(~ind_stayer_mult_cluster))';
    N_clusters = N_clusters + sum(~ind_stayer_mult_cluster);

    % Rename clusters of stayers with multiple clusters
    aux = findgroups(cluster(NT_mover + 1: NT_mover + NT_stayer_mult_cluster), id(NT_mover + 1: NT_mover + NT_stayer_mult_cluster));
    cluster(NT_mover + 1: NT_mover + NT_stayer_mult_cluster) = aux +  N_clusters;
    clear aux

    % Rename cluster ids yet again
    [u_cluster, ~, cluster] = unique(cluster, 'stable');
    N_clusters = numel(u_cluster);
    clear u_cluster
    
    % Generate length of each cluster
    N_C = accumarray(cluster, 1);
    N_C = N_C(cluster); % Bring back to correct length
    
    % Sort clusters by length for movers and stayers separately
    % Sort the data first by cluster size (descending), then by cluster ID (ascending)
    [~, i_sort_size_C_mover] = sortrows([N_C(ind_mover), cluster(ind_mover)], [-1, 2]);    
    [~, i_sort_size_C_stayer] = sortrows([N_C(~ind_mover), cluster(~ind_mover)], [-1, 2]);
    i_sort_size_C = [i_sort_size_C_mover; i_sort_size_C_stayer + length(i_sort_size_C_mover)];

    % Order the rest of the variables
    y = y(i_sort_size_C); id = id(i_sort_size_C); year = year(i_sort_size_C);
    firmid = firmid(i_sort_size_C); cluster = cluster(i_sort_size_C);
    N_C = N_C(i_sort_size_C);
    group = group(i_sort_size_C, :); group_old = group_old(i_sort_size_C, :);
    clear i_sort_size_C  
    
    % Rename identifiers one last time
    [~, ~, firmid] = unique(firmid, 'stable');
    [~, ~, id] = unique(id, 'stable');
    [~, ~, cluster] = unique(cluster, 'stable');


    %%% Count number of workers per firm*year
    firmid_y = findgroups(firmid, year);
    [GC, GR] = groupcounts([id, firmid_y]);

    dict_id_firmid_weight = [GR{:,1} GR{:,2} GC]; % unique id, unique firmid*year and counts per match*year
    % sort per match*year to repeat later on
    [~, i_sort] = sort(dict_id_firmid_weight(:,2)); % sort per firmid_y for repetition
    dict_id_firmid_weight = dict_id_firmid_weight(i_sort,:);

    % Note: Previously need to reorder all the variables (year included) by
    % clusters above
    % Number of observations per firmid*year to decide if LdM moment can be
    % computed and adjust
    N_obs_per_firmid_y = accumarray(dict_id_firmid_weight(:,2),dict_id_firmid_weight(:,3)); % counts of observations per firmid*year
    N_obs_per_firmid_y = N_obs_per_firmid_y(firmid_y); % repeat vector

    clear GR GC dict_id_firmid_weight year i_sort
    if ind_light ==1
        % Clear memory from base workspace
        evalin('base', 'year = [];'); evalin('base', 'year = [];')
    end
    
    % Create sparse matrix of firmid_y to compute sum of worker fixed effects per firmid*year by matrix
    % multiplication
    F_y   = sparse(firmid_y', 1:NT, true); % already the transpose

    % Initialize cell to store unique values of group old (each group may have
    % different size)
    u_group_old = cell(N_multi_group,1);    

    for i=1:N_multi_group
         [~, ~, n] = unique(group(:, i), 'stable');
         group(:, i) = n;
         u_group_old{i} = unique(group_old(:, i), 'stable');
    end
    clear group_old

    % Use cellfun to get the length of each matrix in the cell array
    N_u_group = cellfun(@length, u_group_old);
   

    %% Sort per group for faster covariance computation

    [ch_group, i_group, N_id_group, N_firmid_group, N_obs_group] = compute_group_changes(group, id, firmid, dim_fe, N_multi_group, N_u_group);

    
    %% Compute variance of outcome variable
    %%% Overall and per group variance of log wages
	var_y_overall = var(y);
    % Stack overall variance and variance per group
    % Initialize cell for var_y
    var_y = cell(N_multi_group,1);

    for i =1:N_multi_group
        var_y{i} = [var_y_overall; splitapply(@var,y,group(:,i))];
    end
    clear group


    %% Initial AKM with residualized data. Regression only with movers
    
    % Design matrix first two fixed effects
    D   = sparse(1:NT_mover, id(ind_mover)', 1);
    F   = sparse(1:NT_mover, firmid(ind_mover)', 1);
	
    % Need to normalize one firm effect effect to avoid collinearity
    % Take out one firm effect. Last firm
    F(:, end) = [];

    fprintf('\n')
    disp('---------- Initial AKM ---------- ')    
    X   = [D, -F]; % Only has movers
    xx  = X'*X; % k_mover * k_mover
    xy  = X'*y(ind_mover);
    clear D F

    disp('Building preconditioner for Laplacian Matrix...')
    [~,L] = evalc('cmg_sdd(xx)'); %preconditioner for Laplacian matrices.
    disp('Done!')
        
    % Initial AKM
    b   = pcg(xx, xy, mytol, 300, L); % b only contains worker fixed effect estimates for movers and the firm fixed effect estimates
    clear xy

    %%% ANALYZE RESULTS %%%

    % Estimate worker fixed effects for stayers to compute residuals.
    % Stack vector of fixed effects for movers and stayers
    [pe_tilde, fe_tilde] = compute_fe(b, dim_fe_mover, y, id, firmid, NT_mover);

    % Residuals for the full sample (stayers and movers)
    r   = y - pe_tilde - fe_tilde;     

    %%% Average worker fixed effects at the firm per year
    sum_pe = F_y * pe_tilde;    
    sum_pe = sum_pe(firmid_y); % repeat vector
    
    % Average of estimated worker fixed effects for coworkers
    av_coworkers = (sum_pe - pe_tilde)./(N_obs_per_firmid_y-1); % assign average coworker FE by removing the current worker fixed effect for firms with more than 2 workers
    

    % Plugin estimates in a cell
    plugin = cell(N_multi_group, 1);
    for i= 1: N_multi_group
        pe_tilde = pe_tilde(i_group(:, i));
        fe_tilde = fe_tilde(i_group(:, i));
        av_coworkers = av_coworkers(i_group(:, i));
        plugin{i} = compute_moments_group_LdM(pe_tilde, fe_tilde, ch_group{i}, av_coworkers);
    end

    % Clean memory
    clear pe_tilde fe_tilde b av_coworkers sum_pe

    % Sort results by old group identifier
    u_group_old_sorted = cell(N_multi_group,1);
    i_old_group = cell(N_multi_group,1);

    % Order variables for output
    for i= 1:N_multi_group
        [u_group_old_sorted{i}, i_old_group_aux] = sort(u_group_old{i});
        i_old_group{i} = [1; i_old_group_aux + 1];
        N_obs_group{i} = N_obs_group{i}(i_old_group{i}); % Reorder per ordered old group identifiers
        N_id_group{i} = N_id_group{i}(i_old_group{i});
        N_firmid_group{i} = N_firmid_group{i}(i_old_group{i});

        var_y{i} = var_y{i}(i_old_group{i});
        plugin{i} = plugin{i}(i_old_group{i},:);
    end
	clear i_old_group_aux u_group_old
        
    % Determine label of second order moments
    str_disp = [{'var_worker'}; {'var_firm'}; {'var_av_worker'}; {'cov_worker_firm'}; {'cov_worker_av_worker'}; {'cov_av_worker_firm'};...
        {'N_obs'}; {'N_id'}; {'N_firmid'}];
    keep_no_LdM = [1:2, 4];

    %%% Show in command and, if ind_export=1, write plugin, corrected and N_obs_group files
    write_est_multi_group(ind_export, plugin, 'plugin', str_disp, u_group_old_sorted,...
        N_multi_group, N_obs_group, N_id_group, N_firmid_group, filename, v_filename_group)
    write_N_obs_multi_group(ind_export, u_group_old_sorted,...
        N_multi_group, N_obs_group, N_id_group, N_firmid_group, filename, v_filename_group)

   
    %% VARIANCE-COVARIANCE MATRIX ESTIMATION CORRECTIONS
    fprintf('\n')
    disp('---------- Covariance matrix estimation ---------- ')

    % Pass variables to workers
    X_par  = parallel.pool.Constant(X);
    xx_par = parallel.pool.Constant(xx);
    L_par  = parallel.pool.Constant(@() L);
    
    
    %% Estimation of projection and residual matrices for each cluster
    fprintf('\n')
    disp('---------- Estimation P_gg and M_gg for movers ---------- ')
        
    % Output is the estimated leverage (Pii) and the estimated 1 - leverage (Mii)
    disp('Start to estimate P_gg and M_gg...')
    tic    
    [P_hat, M_hat] = est_M_P(n_lev, NT_mover, X, xx, L, cluster(ind_mover), mytol);
    disp('Done!')
    toc
    %% Estimate leave-cluster-out error for movers

    % For movers with multiple clusters
    disp('Getting leave-cluster-out errors for movers with multiple clusters...')
    tic
    eps_LC = get_eps_LC_cluster(P_hat(1:NT_mover_mult_cluster, 1:NT_mover_mult_cluster),...
        M_hat(1:NT_mover_mult_cluster, 1:NT_mover_mult_cluster), r(1:NT_mover_mult_cluster));
    toc

    % For movers with only one cluster (Reminder: we are using
    % leave-obs out for these workers)
    P_ii = diag(P_hat(NT_mover_mult_cluster + 1 : NT_mover, NT_mover_mult_cluster + 1 : NT_mover));    
    clear P_hat
    M_ii = diag(M_hat(NT_mover_mult_cluster + 1 : NT_mover, NT_mover_mult_cluster + 1 : NT_mover));
    clear M_hat

    inv_Mii_mover_one_cluster = full((P_ii + M_ii)./M_ii);
    clear P_ii M_ii

    eps_LC = [eps_LC; inv_Mii_mover_one_cluster .* r(NT_mover_mult_cluster + 1 : NT_mover)];    
    clear inv_Mii_mover_one_cluster
    
    %% Estimate leave-cluster-out error for stayers

    eps_LC_stayers = get_eps_LC_cluster_stayers(id(~ind_mover), cluster(~ind_mover), N_C(~ind_mover), r(~ind_mover));

    %% Calculate B+ and B-
    
    disp('Calculating B+ and B-...')
    % Start with stayers
    
    disp('for stayers...')
    tic
    % Collapse observations of stayers at the match level: need compute
    % variance of the average errors per stayer
    est_var_r_stayers = collapse_var_stayers_cluster(eps_LC_stayers, y(~ind_mover), dim_fe_mover, ...
        id(~ind_mover), cluster(~ind_mover), N_C(~ind_mover));
    clear eps_LC_stayers

    filt_plus = (est_var_r_stayers>=0);
    % B+
    est_eigenvalue = est_var_r_stayers;
    est_eigenvalue(~filt_plus) = 0;
    B_plus_stayers  = sqrt(est_eigenvalue);
        
    % B-
    est_eigenvalue = abs(est_var_r_stayers);
    est_eigenvalue(filt_plus) = 0;
    B_minus_stayers  = sqrt(est_eigenvalue);
    clear est_eigenvalue est_var_r_stayers
    toc

    disp('and now for movers...')
    tic
    % Now with movers
    [B_plus, B_minus] = decomp_V_LC(eps_LC, y(ind_mover), cluster(ind_mover), N_C(ind_mover));
    toc

    % Stack all the vectors together
    B_plus  = [B_plus; B_plus_stayers]; clear B_plus_stayers;
    B_minus = [B_minus; B_minus_stayers]; clear B_minus_stayers;

     % pass variable to workers
    B_plus_par = parallel.pool.Constant(B_plus);
    B_minus_par = parallel.pool.Constant(B_minus);
    clear B_plus B_minus filt_plus  

    %%   BOOTSTRAP CORRECTION 
    fprintf('\n')
    disp('---------- Bootstrap ---------- ')
    
    %%% Prep objects to do bootstrap and clear unnecessary large objects from memory
    
	% Create vector of firmids for stayers for each id
    id_stayer_norm = id(~ind_mover) - dim_fe_mover(1); % ids of stayers starting from 1
    firmid_stayers = accumarray(id_stayer_norm, firmid(~ind_mover), [], @(x) x(1)); % One firmid per stayer
    cluster_par = parallel.pool.Constant(cluster(ind_mover));
    clear id_stayer_norm ind_mover cluster

    i_group_par = parallel.pool.Constant(i_group);
    clear i_group
    ch_group_par = parallel.pool.Constant(ch_group);
    N_u_group_par = parallel.pool.Constant(N_u_group);
    clear ch_group   

    % pass variable to workers; clear old variables 
    firmid_stayers_par = parallel.pool.Constant(firmid_stayers);
    clear firmid_stayers   
    

    % Vectors of id's and firmid's to assign after having fixed effects
    % estimates
    id_par = parallel.pool.Constant(id);
    clear id
    firmid_par = parallel.pool.Constant(firmid);
    clear firmid
    firmid_y_par = parallel.pool.Constant(firmid_y);
    F_y_par = parallel.pool.Constant(F_y);
    N_obs_per_firmid_y_par = parallel.pool.Constant(N_obs_per_firmid_y);  
    clear firmid ind_mover N_obs_per_firmid_y firmid_y F_y
   
    
    % Run bootstrap
    tic
    delta = bootstrap_cluster_multi_group_LdM(n_boot, B_plus_par.Value, B_minus_par.Value,...
        dim_fe_mover,X_par.Value, xx_par.Value, L_par.Value, mytol, NT_mover,...
        id_par.Value, firmid_par.Value, firmid_stayers_par.Value, cluster_par.Value,...
        i_group_par.Value, ch_group_par.Value,...
        N_u_group_par.Value,...
        N_obs_per_firmid_y_par.Value, firmid_y_par.Value, F_y_par.Value);     
    toc

        %% Write corrected moments and decompositions
    plugin_no_LdM = cell(N_multi_group,1);
    corrected_no_LdM = cell(N_multi_group,1);
    corrected = cell(N_multi_group,1);

    % Order variables for output
    for i = 1:N_multi_group
       delta{i} = delta{i}(i_old_group{i},:);
       corrected{i} = plugin{i} - delta{i};

       plugin_no_LdM{i} = plugin{i}(:, keep_no_LdM);
       corrected_no_LdM{i} = corrected{i}(:, keep_no_LdM);
    end
        
    %%% Show in command and, if ind_export=1, write corrected
    write_est_multi_group(ind_export, corrected, 'corrected', str_disp, u_group_old_sorted,...
        N_multi_group, N_obs_group, N_id_group, N_firmid_group, filename, v_filename_group)
     
    %% Variance decomposition
    fprintf('\n')
    disp('---------- Variance decomposition ---------- ')        
    
    % Determine label of second order moments for the decomposition
    str_disp = [{'var_worker'}; {'var_firm'}; {'2*cov_worker_firm'}]; 

    [decomp_pi, decomp_b] = write_decomp_multi_group(ind_export, plugin_no_LdM, corrected_no_LdM, str_disp, var_y,...
        u_group_old_sorted, dim_fe, N_multi_group, filename, v_filename_group);


    disp('FINISHED!')
end