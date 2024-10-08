function [plugin, delta, corrected, decomp_pi, decomp_b, dim_fe, NT] = correction_LdM(y, id, firmid, year, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                type_hc, type_leave, ind_light, mytol, ind_export, filename)
    %%%% STANDARD CORRECTION %%%%%
    % Description: It only does the bias corrections on the entire sample.
    % It computes the covariance between worker and average co-workers and the variance of average co-worker fixed effects 
    % Suited to handle any diagonal covariance matrix estimate (hom, hc0,
    % hc1, hc2, and hcu).
    % For covariance matrix estimates with clusters at match level use
    % "correction_match.m".
    % For corrections for different groups of the sample use "correction_groups_LdM.m"
    
    % ----------------------------------------------------------------------
    %   START
    % ----------------------------------------------------------------------
    
    %% Define default parameters
    default_param.type_leave = 'obs';
    
    % Assign default values if the variable is missing
    if size(type_leave,2)==0
        type_leave = default_param.type_leave;
    end

    disp(strcat('Setup: leave option is',{' '},type_leave,'. Covariance matrix is assumed diagonal and estimated with ',{' '},type_hc))

    %% Initial variables
    disp(['Size initial data: ',num2str(length(id))])
    
    % Number of FE
    n_FE = 2 + size(other_fe, 2);

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
    y=y(sel); firmid=firmid(sel); id=id(sel); year=year(sel);

    if ~isempty(mat_controls)
        mat_controls = mat_controls(sel,:);  
        if ind_light ==1
            % Clear memory from base and caller workspace
            evalin('base', 'mat_controls = [];'); evalin('base', 'mat_controls = [];')
        end
    end
    
    if n_FE > 2
        other_fe = other_fe(sel,:);
    end

    % Relabel ids after keeping the largest connected set
    disp('Relabeling ids again...')
           
    % Relabel the firms
    [~, ~, n] = unique(firmid);
    firmid = n;

    % Relabel the workers
    [~, ~, n] = unique(id);
    id = n;
    
    % Relabel other Fixed Effects
    if n_FE > 2
        for j = 1 : n_FE - 2
            [~, ~ ,n] = unique( other_fe(:, j) );
            other_fe(:, j) = n;
        end
        if ind_light ==1
            % Clear memory from base workspace
            evalin('base', 'other_fe = [];'); evalin('base', 'other_fe = [];')
        end
    end

    if ind_light ==1
        % Clear memory from base workspace
        evalin('base', 'id = [];'); evalin('base', 'id = [];')
        evalin('base', 'firmid = [];'); evalin('base', 'firmid = [];')
        evalin('base', 'y = [];'); evalin('base', 'y = [];')
    end

    disp(['Size connected set data: ',num2str(size(y,1))])
        
    %% Pruning of data and removal of one-time observations
    
    disp(['Type of data selection procedure: ',type_leave])

    % Indicator of whether pruning is necessary given choice of type_leave 
    ind_pruning = (strcmp(type_leave, 'obs') || strcmp(type_leave, 'match') || strcmp(type_leave, 'worker'));
    
    % Check if chice of variance estimate is incompatible with choice of
    % type_leave
    if (  ( strcmp( type_hc, 'hc_2') || strcmp( type_hc, 'hc_u' ) ) && ind_pruning == 0)
           disp('Current choice of variance estimate (HC_U or HC_2) incompatible with data selection.')
           disp('Proceed using leave-one-out sample')
           type_leave = 'obs';
           ind_pruning = 1;
    end

    if ind_pruning == 1
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
        y = y(sel); firmid = firmid(sel); id = id(sel); year=year(sel);
        if ~isempty(mat_controls)
            mat_controls = mat_controls(sel, :);   
        end        
        if n_FE > 2
            other_fe = other_fe(sel, :);
        end

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
        % Leave observation-out is the default

        if strcmp(type_leave, 'worker') 
            sel = pruning_worker_LdM(id, firmid, year); % leave worker as KSS
        elseif strcmp(type_leave, 'match') 
            sel = pruning_match_LdM(id, firmid, year); % leave match out
        else  
            sel = pruning_obs_LdM(id, firmid, year); % leave observation out
            % sel = pruning_obs_LdM_test(id, firmid, year); % leave observation out
            % sel = pruning_obs_LdM_old(id, firmid, year); % leave observation out
        end

        % Filter variables and rename identifiers again
        y = y(sel); id = id(sel); firmid = firmid(sel);  year=year(sel);
        if ~isempty(mat_controls)
            mat_controls = mat_controls(sel, :);   
        end
        if n_FE > 2
            other_fe = other_fe(sel, :);
        end

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
    end
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

    %% Define stayers and movers to work with them separately
    % Sort movers at the beginning of the sample
    [ind_mover, i_sortmover] = sort_mover(id, firmid); % ind_mover is already sorted (movers come first)
            
    % Order the rest of the variables
    y = y(i_sortmover); id = id(i_sortmover); year = year(i_sortmover);
    firmid = firmid(i_sortmover);
    clear i_sortmover  

    % Rename identifiers one last time
    [~, ~, firmid] = unique(firmid, 'stable');
    [~, ~, id] = unique(id, 'stable');
    
    % Get how many movers are
    NT_mover    = sum(ind_mover);
    dim_fe_mover    = [max(id(ind_mover)) max(firmid(ind_mover))-1];
    disp(['Share of mover ids: ', num2str(dim_fe_mover(1)/dim_fe(1))])

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

    %%% Count number of workers per firm*year
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
    clear GR GC dict_id_firmid_weight year i_sort
    if ind_light ==1
        % Clear memory from base workspace
        evalin('base', 'year = [];'); evalin('base', 'year = [];')
    end

    %%% Get unique ids and firmid
    N_id = numel(unique(id));
    N_firmid = numel(unique(firmid));
    NT = numel(id);

    % Create sparse matrix of firmid_y to compute sum of worker fixed effects per firmid*year by matrix
    % multiplication
    F_y   = sparse(firmid_y', 1:NT, true); % already the transpose

    %%% Average worker fixed effects at the firm per year
    sum_pe = F_y * pe_tilde;    
    sum_pe = sum_pe(firmid_y); % repeat vector
    
    % Average of estimated worker fixed effects for coworkers
    av_coworkers = (sum_pe - pe_tilde)./(N_obs_per_firmid_y-1); % assign average coworker FE by removing the current worker fixed effect for firms with more than 2 workers

 
    % Plugin estimates including LdM moment
    plugin = compute_moments_LdM(pe_tilde, fe_tilde, av_coworkers);

    % Clean memory
    clear pe_tilde fe_tilde b sum_pe  av_coworkers

    % Determine label of second order moments for the decomposition
    str_disp = [{'var_worker'}; {'var_firm'}; {'var_av_worker'}; {'cov_worker_firm'}; {'cov_worker_av_worker'}; {'cov_av_worker_firm'}; {'N_obs'}; {'N_id'}; {'N_firmid'}];
    pos_LdM = [3,5:6]; 
    
    % Write results
    write_est(ind_export, plugin, 'plugin', str_disp, NT, N_id, N_firmid, filename)
   
    %% VARIANCE-COVARIANCE MATRIX ESTIMATION CORRECTIONS
    fprintf('\n')
    disp('---------- Covariance matrix estimation ---------- ')

    % Pass variables to workers
    X_par  = parallel.pool.Constant(X);
    xx_par = parallel.pool.Constant(xx);
    L_par  = parallel.pool.Constant(@() L);
    
    % Clean some memory
    clear xx X L

    if strcmp(type_hc, 'hom')
		est_var_r = (1 ./ dof) .* sum(r.^2) .* ones(NT,1);
    elseif strcmp(type_hc, 'hc_0')
        est_var_r = r.^2;
    elseif strcmp(type_hc, 'hc_1')
        hc = NT ./ dof;
        est_var_r = r.^2.*hc;
    elseif ( strcmp( type_hc, 'hc_2') || strcmp( type_hc, 'hc_u' ) )
        
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

        if (strcmp(type_hc,'hc_2'))
		    est_var_r = r.^2.*inv_Mii; % This is hc_2 in the paper      
            
            % Check if remaining issues
            if (sum(est_var_r<0)>0)
               disp('THERE ARE REMAINING ISSUES') 
            end
        else 
            est_var_r = y.*r.*inv_Mii; % This is hc_u in the paper 
        end   
    else 
        disp('Unrecognized type of variance estimate')
    end
    
    clear r Mii inv_Mii inv_Mii_mover correction_lev

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

    if strcmp(type_hc,'hc_u')
        % B-
        est_eigenvalue = abs(est_var_r);
        est_eigenvalue(filt_plus) = 0;
        
        B_minus  = sqrt(est_eigenvalue);
        % pass variable to workers 
        B_minus_par = parallel.pool.Constant(B_minus);
        clear B_minus est_eigenvalue filt_plus est_var_r
    else
        B_minus = [];
        % pass variable to workers 
        B_minus_par = parallel.pool.Constant(B_minus);
        clear B_minus est_eigenvalue filt_plus est_var_r
    end

    %%   BOOTSTRAP CORRECTION 
    fprintf('\n')
    disp('---------- Bootstrap ---------- ')
    
    %%% Prep objects to do bootstrap and clear unnecessary large objects from memory
    % Variance of outcome variable
    var_y =var(y); % variance of log wages removing the observations where we cannot compute the LdM moment
    clear y 
    
	% Create vector of firmids for stayers for each id
    id_stayer_norm = id(~ind_mover) - dim_fe_mover(1); % ids of stayers starting from 1
    firmid_stayers = accumarray(id_stayer_norm, firmid(~ind_mover), [], @(x) x(1)); % One firmid per stayer
    clear id_stayer_norm ind_mover

    % pass variable to workers; clear old variables 
    firmid_stayers_par = parallel.pool.Constant(firmid_stayers);
    clear firmid_stayers   
   
    % Vectors of id's and firmid's to assign after having fixed effects
    % estimates
    id_par = parallel.pool.Constant(id);
    clear id
    firmid_par = parallel.pool.Constant(firmid);
    firmid_y_par = parallel.pool.Constant(firmid_y);
    F_y_par = parallel.pool.Constant(F_y);
    N_obs_per_firmid_y_par = parallel.pool.Constant(N_obs_per_firmid_y);  
    clear N_obs_per_firmid_y firmid firmid_y F_y
    
    % Run bootstrap
    tic
    [delta, corrected] = bootstrap_LdM(n_boot, B_plus_par.Value, B_minus_par.Value,...
        dim_fe_mover,X_par.Value, xx_par.Value, L_par.Value, plugin, mytol, NT_mover,...
        id_par.Value, firmid_par.Value, firmid_stayers_par.Value,...
        N_obs_per_firmid_y_par.Value, firmid_y_par.Value, F_y_par.Value);     
    toc

    % Write results
    write_est(ind_export, corrected, 'corrected', str_disp, NT, N_id, N_firmid, filename)
    

    %% Variance decomposition
    fprintf('\n')
    disp('---------- Variance decomposition ---------- ')
         
    %%% Estimates without LdM moment 
    plugin_no_LdM = plugin;
    plugin_no_LdM(:,pos_LdM) = []; % remove moments related to LdM
    corrected_no_LdM = corrected;
    corrected_no_LdM(:,pos_LdM) = []; % remove moments related to LdM

    str_disp = [{'var_worker'}; {'var_firm'}; {'2*cov_worker_firm'}];
    
    [decomp_pi, decomp_b] = write_decomp(ind_export, plugin_no_LdM, corrected_no_LdM, str_disp, var_y,...
        dim_fe, filename);


    disp('FINISHED!')
end