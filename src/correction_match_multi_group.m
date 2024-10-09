function [plugin, delta, corrected, decomp_pi, decomp_b, dim_fe, NT] = correction_match_multi_group(y, id, firmid, group, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                type_leave, ind_light, mytol, ind_export, filename, v_filename_group)
    %%%% CORRECTION FOR DIFFERENT GROUPS WITH LEAVE-MATCH-OUT VARIANCE ESTIMATE %%%%%
    % Description: This code does the bias corrections on the entire sample and
    % different subsamples as determined by the group vector.
    % Suited for covariance matrix estimate with clusters at the match level
    % For diagonal covariance matrix estimate (hom, hc0, hc1, hc2, and hcu)
    % use "correction_group.m"
    
    % ----------------------------------------------------------------------
    %   START
    % ----------------------------------------------------------------------
    
    %% Define default parameters
    default_param.type_leave = 'match';
        
    % Assign default values if the variable is missing
    if size(type_leave,2)==0
        type_leave = default_param.type_leave;
    end

    disp(strcat('Setup: leave option is',{' '},type_leave,'. Covariance matrix is assumed block diagonal'))

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
    [~, ~, n] = unique(id);
    id = n;

    % Relabel other Fixed Effects
    if n_FE > 2
        for j = 1 : n_FE - 2
            [~, ~ ,n] = unique( other_fe(:,j) );
            other_fe(:,j) = n;
        end
    end
            

	%% Compute the connected set and relabel
    disp('---------- Connected set ---------- ')

    % Find the connected set. Then define dimensions and relevant matrices.
    sel = connected(firmid, lagfirmid, sel);

    % Filter after finding connected set
    y=y(sel); firmid=firmid(sel); id=id(sel);

    if ~isempty(mat_controls)
        mat_controls = mat_controls(sel,:);  
        if ind_light ==1
            % Clear memory from base and caller workspace
            evalin('base', 'mat_controls = [];'); evalin('caller', 'mat_controls = [];')
        end

    end    
    if n_FE > 2
        other_fe = other_fe(sel,:);
        if ind_light ==1
            % Clear memory from base and caller workspace
            evalin('base', 'other_fe = [];'); evalin('caller', 'other_fe = [];')
        end

    end
    group = group(sel,:); 
    group_old = group_old(sel,:);

    % Relabel the workers
    [~, ~, id] = unique(id);
       
    N_multi_group = size(group,2);

    % Update group
    for i=1:N_multi_group
         [~, ~, n] = unique(group(:, i));
         group(:, i) = n;
    end

    if ind_light ==1
        % Clear memory from base and caller workspace
        evalin('base', 'id = [];'); evalin('caller', 'id = [];')
        evalin('base', 'firmid = [];'); evalin('caller', 'firmid = [];')
        evalin('base', 'group = [];'); evalin('caller', 'group = [];')
        evalin('base', 'y = [];'); evalin('caller', 'y = [];')
    end


    disp(['Size connected set data: ',num2str(size(y,1))])
        
    %% Pruning of data and removal of one-time observations
    
    disp(['Type of data selection procedure: ', type_leave])

    % Indicator of whether the right choice of pruning was chosen 
    ind_pruning = (strcmp(type_leave, 'match') || strcmp(type_leave, 'worker'));
    
    if ind_pruning == 0
       disp('Current choice of variance estimate (Leave-match-out) incompatible with chosen data selection method.')
       disp('Proceed using leave-match-out sample') 
       type_leave = 'match';
    end

    disp('---------- Pruning ---------- ')   
    tic           

    % Remove ids that appear only once. They dont generate connections
    % and have leverage equal to 1

    % Count obs per id
    [~, ~, iunique] = unique(id);
    n = accumarray(iunique(:), 1);

    % Remove observations with only one observation
    sel = n > 1; 
    sel = sel(iunique); % Repeat to filter original vector
    clear iunique n
           
    % Filter out one-time ids
    y = y(sel); firmid = firmid(sel); id = id(sel);
    if ~isempty(mat_controls)
        mat_controls = mat_controls(sel, :);   
    end        
    if n_FE > 2
        other_fe = other_fe(sel, :);
    end
    group = group(sel,:);
    group_old = group_old(sel,:);

    disp(['Size after removing one-timers: ', num2str(size(y,1))])
              
    % Rename FE ids after removing one-timers
    [~, ~, firmid] = unique(firmid, 'stable');
    [~, ~, id] = unique(id, 'stable');
    if n_FE > 2
       for j = 1 : n_FE - 2
           [~, ~ ,n] = unique( other_fe(:, j), 'stable' );
           other_fe(:, j) = n;
       end
    end

    % Choose pruning type: leave worker, match or observation out
    % Leave match is the default

    if strcmp(type_leave, 'worker') 
        sel = pruning_worker(id, firmid); % leave worker as KSS
    else  
        sel = pruning_match(id, firmid); % leave match out
    end

    % Filter variables and rename identifiers again
    y = y(sel); id = id(sel); firmid = firmid(sel);
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
    if n_FE > 2
       for j = 1 : n_FE - 2
           [~, ~ ,n] = unique( other_fe(:, j), 'stable' );
           other_fe(:, j) = n;
       end
    end

    toc
        
    disp(['Size of pruned data: ', num2str(size(y,1))])
    
    clear lagfirmid sel n    
    
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
    
    %% Define stayers and movers to work with them separately. Define clusters at the match level
    % Sort data. Movers first, in descending order by observations per
    % cluster. Then stayers.        
    [ind_mover, cluster, i_sortmover, ~] = sort_mover_match(id, firmid); % ind_mover is already sorted
       
    % Order the rest of the variables
    y = y(i_sortmover); id = id(i_sortmover);
    firmid = firmid(i_sortmover);
    group = group(i_sortmover,:); group_old = group_old(i_sortmover,:);
    clear i_sortmover  

    % Rename identifiers one last time
    [~, ~, firmid] = unique(firmid, 'stable');
    [~, ~, id] = unique(id, 'stable');
    
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
    
    % Get how many movers are
    NT_mover    = sum(ind_mover);
    dim_fe_mover    = [max(id(ind_mover)) max(firmid(ind_mover))-1];
    disp(['Share of mover ids: ', num2str(dim_fe_mover(1)/dim_fe(1))])
    
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

    %% Collapse observations (for movers) at the cluster level
    disp('Collapsing observations of movers at the cluster level')
    tic
    id_old = id;
    firmid_old = firmid;
    [y, id, firmid, weights, ind_mover] = collapse_cluster(y, id, firmid, cluster, dim_fe_mover);
    toc

    %% Renormaize dimensions after collapse at cluster level
    NT_original = NT;
    NT_mover_original = NT_mover;

    NT = length(id);
    NT_mover    = sum(ind_mover);

    % Design matrix first two fixed effects
    D   = sparse(1:NT_mover, id(ind_mover)', 1);
    F   = sparse(1:NT_mover, firmid(ind_mover)', 1);
	
    % Need to normalize one firm effect effect to avoid collinearity
    % Take out one firm effect. Last firm
    F(:, end) = [];

    %% Initial AKM with residualized data. Regression only with movers
    fprintf('\n')
    disp('---------- Initial AKM ---------- ')    
    X  = [D, -F]; % Only has movers
    xy = X' * (weights(ind_mover) .* y(ind_mover)); % This is the same as doing X*W*y, where W is a diagonal weighting matrix
    X  = ( (weights(ind_mover)).^0.5 ) .* X; % Add weights to solve the correct normal equations  
    xx = X' * X; % k_mover * k_mover
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
    
    % Weight outcome variable to compute residuals and variance estimates
    y = ( (weights).^0.5 ).*y; % Stayers have weight of 1

    % Residuals for the full sample (stayers and movers)
    r = y - ( (weights).^0.5 ).*(pe_tilde + fe_tilde);     
    
    %%% Plugin estimates

    % Recompute estimated fixed effects, repeating them to have the same
    % observations as the original sample before collapsing at the cluster
    % level
    [pe_tilde, fe_tilde] = compute_fe(b, dim_fe_mover, y, id_old, firmid_old, NT_mover);
    % Plugin estimates in a cell
    plugin = cell(N_multi_group, 1);
    for i= 1: N_multi_group
        pe_tilde = pe_tilde(i_group(:, i));
        fe_tilde = fe_tilde(i_group(:, i));
        plugin{i} = compute_moments_group(pe_tilde, fe_tilde, ch_group{i});
    end

    % Clean memory
    clear pe_tilde fe_tilde b weights
    
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
    str_disp    = [{'var_worker'}; {'var_firm'}; {'cov_worker_firm'}; {'N_obs'}; {'N_id'}; {'N_firmid'}];

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
    
    % Clean some memory
    clear xx X L
    
    %% Estimation of leverage
    fprintf('\n')
    disp('---------- Leverage estimation ---------- ')
        
    % Output is the estimated leverage (Pii) and the estimated 1 - leverage (Mii)
    tic
    disp('Start to estimate the leverages...')
    [~, Mii, correction_lev] = est_lev(n_lev, NT_mover, X_par.Value, xx_par.Value, L_par.Value, id, mytol);
    disp('Done!')
    toc
        
    %% Diagnostic of leverage estimation for movers
    fprintf('\n')
    disp('---------- Diagnostic ---------- ')
    tic
    [inv_Mii_mover, ~] = diagnostic_lev(correction_lev, NT_mover, X_par.Value, xx_par.Value, L_par.Value, Mii(ind_mover), mytol);
    toc
        
    inv_Mii = [inv_Mii_mover; 1./Mii(~ind_mover)];

    est_var_r = y.*r.*inv_Mii; % This is hc_u in the paper
    
    clear r Mii inv_Mii y

    %% Calculate B+ and B-
    
    % Collapse observations of stayers at the match level: need compute
    % variance of the average errors per stayer
    [est_var_r] = collapse_var_stayers(est_var_r, dim_fe_mover, id);

    % Covariance matrix is diagonal. B is just sqrt of estimated
    % variance of the residual because Q is the identity.
    filt_plus = (est_var_r>=0);

    % B+
    est_eigenvalue = est_var_r;
    est_eigenvalue(~filt_plus) = 0;   
    
    B_plus  = sqrt(est_eigenvalue);
    % pass variable to workers
    B_plus_par = parallel.pool.Constant(B_plus);
    clear B_plus est_eigenvalue

    % B-
    est_eigenvalue = abs(est_var_r);
    est_eigenvalue(filt_plus) = 0;
    
    B_minus  = sqrt(est_eigenvalue);
    % pass variable to workers 
    B_minus_par = parallel.pool.Constant(B_minus);
    clear B_minus est_eigenvalue est_var_r filt_plus      

    %%   BOOTSTRAP CORRECTION 
    fprintf('\n')
    disp('---------- Bootstrap ---------- ')
           
	%%% Prep objects to do bootstrap and clear unnecessary large objects from memory
	% Create vector of firmids for stayers for each id
    id_stayer_norm = id(~ind_mover) - dim_fe_mover(1); % ids of stayers starting from 1
    firmid_stayers = accumarray(id_stayer_norm, firmid(~ind_mover), [], @(x) x(1)); % One firmid per stayer
    clear id_stayer_norm ind_mover id firmid

    % pass variable to workers; clear old variables 
    firmid_stayers_par = parallel.pool.Constant(firmid_stayers);
    clear firmid_stayers   
    i_group_par = parallel.pool.Constant(i_group);
    clear i_group
    ch_group_par = parallel.pool.Constant(ch_group);
    N_u_group_par = parallel.pool.Constant(N_u_group);
    clear ch_group   

    % Important to use id and firmid BEFORE collapse at the cluster level
    % Use these assignments to repeat estimates to sample before collapsing
    % to estimate moments for different groups
    id_par = parallel.pool.Constant(id_old);
    clear id_old
    firmid_par = parallel.pool.Constant(firmid_old);
    clear firmid_old
       
    % Run bootstrap
    tic
    delta = bootstrap_multi_group(n_boot, B_plus_par.Value, B_minus_par.Value,...
        dim_fe_mover, X_par.Value, xx_par.Value, L_par.Value, mytol, NT_mover,...
        id_par.Value, firmid_par.Value, firmid_stayers_par.Value, i_group_par.Value, ch_group_par.Value,...
        N_u_group_par.Value);     
    toc

    %% Write corrected moments and decompositions
    % Order variables for output
    corrected = cell(N_multi_group,1);
    for i= 1:N_multi_group
        delta{i} = delta{i}(i_old_group{i},:);
        corrected{i} = plugin{i} - delta{i};
    end
        
    %%% Show in command and, if ind_export=1, write plugin, corrected and N_obs_group files
    write_est_multi_group(ind_export, corrected, 'corrected', str_disp, u_group_old_sorted,...
        N_multi_group, N_obs_group, N_id_group, N_firmid_group, filename, v_filename_group)
  
    %% Variance decomposition
    fprintf('\n')
    disp('---------- Variance decomposition ---------- ')        
    
    % Determine label of second order moments for the decomposition
    str_disp = [{'var_worker'}; {'var_firm'}; {'2*cov_worker_firm'}]; 

    [decomp_pi, decomp_b] = write_decomp_multi_group(ind_export, plugin, corrected, str_disp, var_y,...
        u_group_old_sorted, dim_fe, N_multi_group, filename, v_filename_group);
    
    disp('FINISHED!')
end