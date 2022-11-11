function [plugin,delta,corrected,decomp_pi,decomp_b,dimensions,NT,n_problems_lev]=boot_correction(y,mat_id_fe,mat_controls,n_lev,n_boot,period,group,type_boot,type_hc,mytol,resid_controls,filename)

%{
  
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
    %% 					GENERAL DESCRIPTION
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-  

    This function computes the bootstrap correction of second order
    moments of multi-way fixed effects or two-way fixed effects (e.g. in the labor market
    application, a decomposition of log wages into worker and firm fixed
    effects). We use AKM jargon (workers, firms) when describing the code for
    simplicity.
    
    The mandatory input is a person-year dataset that has to be sorted
    by workers' identifiers (id) and year. The function requires as input
    the log wages or the outcome variable and a matrix with worker and firm 
    identifiers. The rest of the parameters take default values if they are not provided.

    The function automatically performs computation of the largest connected set and leave
    out connected set and computes corrected second order moments and a
    variance decomposition
    
    % Version:
    1.0: First version. 27/09/2020.
    2.0: Replication package. 04/10/2022.
    2.1: Change preconditioner for regressions with two leading fixed effects.
    
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
    %% 					DESCRIPTION OF THE INPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                            %-MANDATORY INPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

    . y: Outcome. Dimensions: N* x 1; N*= # of person-year observations.

    . mat_id_fe: Matrix with fixed effect identifiers. Requires at least
        having worker and firm identifiers. Dimensions: N* x 2
        
        - id: worker indicators. Dimensions: N* x 1
        - firmid: firm indicators. Dimensions: N* x 1


    
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                        %---NON-MANDATORY INPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
    
    . mat_id_fe: The user can add more fixed effects than worker and firm
    ones but the connected set and the leave-one-out connected set are
    computed using id and firmid. The additional fixed effect identifiers
    are normalized to avoid multicolinearities.

    . mat_controls: Matrix of controls. Dimensions: N* x K. The variances
    and covariances are estimated for all K variables together.
    
    . n_lev: Number of simulations for leverage estimation. Mandatory when
    choosing type_hc="hc_2".  A natural number. The default value is 300.

    . n_boot: Number of bootstrap simulations. This parameter governs the
    precision of the estimation. A natural number. The default value is 300.

    . period: Vector stating if the estimation of the second order moments
    and the variance decomposition needs to be splitted into different
    periods. Dimensions:  N* x 1

    . group: Group identifiers when assuming that errors are serially
    correlated. It requires that type_boot='block'. The group can be the worker-firm match, the region or any
    other grouping. Dimensions:  N* x 1

    . type_boot: Can take values of 'diag' when assuming a diagonal
    covariance matrix or 'wild_block' when assuming serial correlation of the
    error terms. The default value is 'diag'.

    . type_hc: Can take values of 'hom', 'hc_0', 'hc_1' or 'hc_2' depending on the
    estimator of the covariance matrix. The default value is 'hc_2'.

    . mytol: Tolerance for pcg when solving the normal equations. The
    default value is 1e-6.

    . resid_controls: Indicator of whether to residualize all other
    variables except the two leading fixed effects. The variance
    decomposition would be done over that residualized left hand side
    variable. For example consider a model y = id + firmid + X_c*beta_c + u.
    If resid_controls == 1 then the first model is estimated and a new
    y_resid is defined as y_resid = y - X_c*beta_hat_c. The default value
    is to residualize with resid_controls=1.

    . filename: path and name to store the second order moment estimates and
    the variance decompositions. String.
    If the user provides filename = 'example' the output files will be: 
        example_plugin_estimates.csv 
        example_corrected_estimates.csv
        example_var_decomp_plugin.csv
        example_var_decomp_corrected.csv




    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
    %% 					DESCRIPTION OF THE OUTPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
 
    . plugin: Plugin estimates of second order moments

    . delta: Bootstrap estimate of the bias of second order moments

    . corrected: Bootstrap corrected second order moments

    . decomp_pi: Variance decomposition using the plugin estimates of
    second order moments. If the user provided a 'period' vector, the
    decomposition is performed per period.

    . decomp_b: Variance decomposition using the bootstrap corrected
    estimates of second order moments. If the user provided a 'period' vector, the
    decomposition is performed per period.

    . dimensions: Dimensions taking into account the normalizations

    . NT: # of person-year observations in the final sample where the
    second order moments are computed. NT can be below N* if the user did
    not provide a connected set or a leave-one-out connected set.

    . n_problems_lev: Number of problematic leverage estimates identified
    in the diagnostic and directly computed.

    . Stored CSVs: The function stores 4 csv files with the plugin
    estimates, the corrected estimates and the respective variance
    decompositions. The user can choose the starting filename for the
    output with the optional input 'filename'.

%}

    %% Initial variables

    % ----------------------------------------------------------------------
    %   START
    % ----------------------------------------------------------------------
    
    % Define default parameters
    default_param.controls = [];
    default_param.n_lev = 300;
    default_param.n_boot = 300;
    default_param.period = [];
    default_param.group = [];
    default_param.type_boot = 'diag';
    default_param.type_hc = 'hc_2';
    default_param.mytol = 1e-6;
    default_param.resid_controls = 1;
    default_param.filename = '';
    
    % Assign default values if the variable is missing
    if nargin < 2
        error('More arguments needed');
    else
        if size(mat_controls,2)==0
            mat_controls = default_param.controls;
        end
        if size(n_lev,2)==0
            n_lev = default_param.n_lev;
        end
        if size(n_boot,2)==0
            n_boot = default_param.n_boot;
        end
        if size(period,2)==0
            period = default_param.period;
        end            
        if size(group,2)==0
            group = default_param.group;
        end  
        if size(type_boot,2)==0
            type_boot = default_param.type_boot;
        end              
        if size(type_hc,2)==0
            type_hc = default_param.type_hc;
        end           
         if size(mytol,2)==0
            mytol = default_param.mytol;
        end           
        if size(resid_controls,2)==0
            resid_controls = default_param.resid_controls;
        end      
        if size(filename,2)==0
            filename = default_param.filename;
        end      
    end
    
    % Required parameters for some cases
    if (size(group,2)==0 && (strcmp(type_boot,'wild_block')))
        disp('Warning: You did not provide a group identifier for the Wild block bootstrap. The code will proceed assuming serial correlation at the match level')
        group = fidgroups(mat_id_fe(:,1),mat_id_fe(:,2));
    end

    % ----------------------------------------------------------------------
    %   GENERATE VECTORS
    % ----------------------------------------------------------------------
    if exist('period')==0 || isempty(period)
       u_period = 1;
    else
       u_period = unique(period);
    end

    %%% First column of id_fe is the worker id and second one firmid.
    id = mat_id_fe(:,1);
    firmid = mat_id_fe(:,2);
    
    disp(['Size initial data: ',num2str(length(id))])
    
    % Stack all the covariates without collinearities. The first two FE will
    % be labeled explicitly. The rest of FE and other controls
    % filtered/renamed within all_covariates matrix. Importantly the
    % connected set or the leave-one-out connected set will only be
    % estimated with the first two fixed effects.
    all_covariates = [mat_id_fe mat_controls];

    
    % Number of FE
    n_FE = size(mat_id_fe,2);
    
    % Design matrices
    id_old = id;
    firmid_old = firmid;
    
    %Lagfirmid
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker

    % Rename FE identifiers
    disp('Relabeling ids...')
    NT = length(y);
    sel=~isnan(lagfirmid);

    %relabel the firms
    [firms,m,n]=unique([firmid;lagfirmid(sel)]);

    firmid=n(1:NT);
    lagfirmid(sel)=n(NT+1:end);
    
    %relabel the workers
    [ids,m,n]=unique(id);
    id=n;
    
    % relabel other FE
    if n_FE>2
        for j=3:n_FE
            [~,~,n]=unique(all_covariates(:,j));
            all_covariates(:,j)=n;
        end
    end

    disp('---------- Connected set ---------- ')
    
	%% Compute the connected set and relabel
    %Find the connected set. Then define dimensions and relevant matrices.
    sel = connected(firmid,lagfirmid,sel);

    %%% Filter after finding connected set
    y=y(sel); firmid=firmid(sel); id=id(sel); 
    all_covariates = all_covariates(sel,:);    
    lagfirmid=lagfirmid(sel);
    if ~isnan(period)
        period=period(sel); % if there is a grouping variable for the corrections update
    end
    if ~isnan(group)
        group=group(sel); % if there is a grouping variable for the corrections update
    end
    
    id_old=id_old(sel); firmid_old=firmid_old(sel);
    
    % Rename FE ids
    disp('Relabeling ids again...')
    NT = size(y,1);
    sel=~isnan(lagfirmid);

    
    %relabel the firms
    [~,~,n]=unique([firmid;lagfirmid(sel)]);

    firmid=n(1:NT);
    lagfirmid(sel)=n(NT+1:end);

    %relabel the workers
    [~,~,n]=unique(id);
    id=n;
    
    % relabel other FE
    if n_FE>2
        for j=3:n_FE
            [~,~,n]=unique(all_covariates(:,j));
            all_covariates(:,j)=n;
        end
    end
    
    disp(['Size connected set data: ',num2str(size(y,1))])
    
    % update mat_id_fe after filters to correctly compute FE dimensions
    mat_id_fe = all_covariates(:,1:n_FE);
    
    % update period
    if ~isnan(period)
        period = findgroups(period); % if there is a grouping variable for the corrections update
    end
    
    
    
    %% Pruning of data and removal of one-time observations
    if (strcmp(type_boot,'diag') && (strcmp(type_hc,'hc_2')))
        disp('---------- Pruning ---------- ')     
        tic
        [y,id,firmid,all_covariates,period,id_old,firmid_old] = pruning(y,id,firmid,all_covariates,n_FE,period,id_old,firmid_old);
        toc
        
        disp(['Size of pruned data: ',num2str(size(y,1))])
        
        
        % Remove ids that appear only once. They dont generate connections
        % and have leverage equal to 1
        % Count obs per id
        [~,~,iunique] = unique(id);
        n = accumarray(iunique(:),1);
        % Remove observations with only one observation
        sel = (n>1); % filter out unique
        sel = sel(iunique); % Repeat to filter original vector
        
        % Filter out one-time ids
        y=y(sel); firmid=firmid(sel); id=id(sel); 
        all_covariates = all_covariates(sel,:);
        if ~isnan(period)
            period=period(sel); % if there is a grouping variable for the corrections update
        end
        lagfirmid=lagfirmid(sel);
        id_old=id_old(sel); firmid_old=firmid_old(sel);
       
        disp(['Size after removing one-times: ',num2str(size(y,1))])
        
        % Rename FE ids after removing one times
        NT = size(y,1);
        sel=~isnan(lagfirmid);

        %relabel the firms
        [~,~,n]=unique([firmid;lagfirmid(sel)]);

        firmid=n(1:NT);
        lagfirmid(sel)=n(NT+1:end);

        %relabel the workers
        [~,~,n]=unique(id);
        id=n;
        
        % relabel other FE in all_covariates
        if n_FE>2
            for j=3:n_FE
                [~,~,n]=unique(all_covariates(:,j));
                all_covariates(:,j)=n;
            end
        end
        
             
        % update mat_id_fe after pruning to correctly compute FE dimensions
        mat_id_fe = all_covariates(:,1:n_FE);
      
        % update period
        if ~isnan(period)
            period = findgroups(period); % if there is a grouping variable for the corrections update
        end
    
         
    end
    
    
    %% Normalized dimensions, design matrices
    NT = length(id);
    % Raw and normalized dimensions of FE
    raw_dim_fe = arrayfun( @(c) numel( unique( mat_id_fe(:, c) )), 1:size( mat_id_fe, 2 )) ; %Get number of unique identifiers for each fixed effect
    dim_fe = raw_dim_fe;
    dim_fe(2:end) = dim_fe(2:end)-1; % Normalizations (substract 1 identifier for all fixed effects except the first one to avoid collinearity)
    
    %Normalized dimensions
    % Vector with dimensions of parameters. Normalized FE dimensions
    if exist('mat_controls') && isempty(mat_controls)==0
        dim_other = size(mat_controls,2);
        %Normalized dimensions
        dimensions = [dim_fe dim_other];
    else
       %Normalized dimensions
        dimensions = dim_fe;
    end
    
    
    %%% Design matrix first two fixed effects
    D=sparse(1:NT,id',1);
    F=sparse(1:NT,firmid',1);

   % Design matrix other FE in all_covariates
    if n_FE>2
        % cell with design matrices
        design_other=cell(n_FE-2,1);
        for j=3:n_FE
            design_other{j-2} = sparse(1:NT,all_covariates(:,j)',1);
        end
    end
    

    % Need to normalize one firm effect and one worker effect to avoid
    % collinearity of fixed effects (because we have fixed effects per
    % occupation)

    % Normalization matrices
    % Take out one firm effect. Last firm
    S=speye(dimensions(2));
    S=[S;sparse(-zeros(1,dimensions(2)))];  %Jx(J-1) matrix for normalization firm effect

    if n_FE>2
        % cell with design matrices
        norm_other=cell(n_FE-2,1);
        for j=3:n_FE
            norm_other{j-2} = speye(dimensions(j));
            norm_other{j-2} = [norm_other{j-2};sparse(-zeros(1,dimensions(j)))]; %Dimensions has normalized dimensions
        end
        
        % cell with normalized other FE
        cell_other = cellfun(@(x,y) x*y, design_other,norm_other, 'UniformOutput',false);
        X_other_fe = reshape(cell2mat(cell_other),NT,[]);
    end
    
	%% Decide if it is possible to have a Laplacian matrix representation
    %%% Matrix with all the covariates for the regression and cell with
    %%% design matrices
    lap=0;
    if n_FE>2 && exist('mat_controls') && isempty(mat_controls)==0 
        
        design_mat_first = [{D} {F*S}];
        design_mat = [design_mat_first {X_other_fe}  {all_covariates(:,n_FE+1:end)}];
         
        clear design_mat_first
                
        X=[D,F*S,X_other_fe,all_covariates(:,(n_FE+1):end)];       % our FE without multicollinearities and other covariates
        
    elseif n_FE>2 && (exist('mat_controls')==0 || isempty(mat_controls)==1)
        
        design_mat_first = [{D} {F*S}];
        design_mat = [design_mat_first X_other_fe];
        
        clear design_mat_first
        
        X=[D,F*S,X_other_fe];       % our FE without multicollinearities
        
    elseif n_FE==2 && exist('mat_controls')  && isempty(mat_controls)==0 
        
        design_mat_first = [{D} {F*S}];
        design_mat = [design_mat_first {all_covariates(:,(n_FE+1):end)}];
        
        clear design_mat_first
        
        X=[D,F*S,all_covariates(:,(n_FE+1):end)];       % our FE without multicollinearities and other covariates
    else    
        
        lap=1; %indicator that X'X can be written as a Laplacian matrix
    end
    
    %% Residualize extra Fixed Effects and Controls (if resid_controls==1) 
    if lap==1 && resid_controls==1
        disp('Warning: no controls nor extra fixed effects to residualize...')
        disp('Proceeding with two leading fixed effects')
    elseif lap==0 && resid_controls==1
        disp('Residualizing controls and/or extra fixed effects')
        n_2FE = size([D,F*S],2); %dimension of two leading fixed effects
        xx=X'*X;
        xy=X'*y;
        % Preconditioner for general matrices
        L =lchol_iter(xx);
       if size(L,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
           % Estimation with all covariates
                b           = pcg(xx,xy,mytol,300,L,L');
       else 
                b           = pcg(xx,xy,mytol,300);
       end
        % Redefine left-hand variable. Variance decompositions will be done 
        %on this variable
        y = y - X(:,n_2FE+1:end)*b(n_2FE+1:end);

        %Turn on that the remaining covariate matrix can form a Laplacian
        %design matrix.
        lap = 1;
        dimensions = dim_fe(1:2);
        disp('Done!')
    end
    
    %% Initial AKM
    disp('---------- Initial AKM ---------- ')
    
    if lap==1
        X=[D,-F*S];
        design_mat = [{D} {-F*S}];
        xx=X'*X;
        xy=X'*y;
        disp('Building preconditioner for Laplacian Matrix...')
        [results,L] = evalc('cmg_sdd(xx)'); %preconditioner for Laplacian matrices.
        disp('Done!')
            
        % Initial AKM
        b = pcg(xx,xy,mytol,300,L);
    else
        xx=X'*X;
        xy=X'*y;
        disp('Building preconditioner for Non-Laplacian Matrix...')
        L=lchol_iter(xx);
        disp('Done!')
        
        % Initial AKM
        if size(L,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
            % Estimation with all covariates
            b = pcg(xx,xy,mytol,300,L,L');
        else 
            b = pcg(xx,xy,mytol,300);
        end
    end
    
    % ANALYZE RESULTS
    xb0=X*b;     % fitted values
    r=y-xb0;     % residuals

    % Degrees of freedom
    dof=NT-sum(dimensions);  
    RMSE=sqrt(sum(r.^2)/dof);

    %%% Cell with betas. Flexible for many fixed effects and covariates
    betas = cell(1,size(design_mat,2));
    for i=1:size(design_mat,2)
        if i==1
            pos = 0;
        else
            pos = pos + dimensions(i-1);
        end
        betas{i} = b((pos+1):(pos+dimensions(i)));
    end 
    
    cell_all_fitted = cellfun(@(x,y) x*y, design_mat,betas, 'UniformOutput',false);
    
    % Convert cell to matrix with fitted values in columns
    all_fitted = reshape(cell2mat(cell_all_fitted),NT,[]);
    
    
    %%% Filter for off-diagonal elements
    if length(u_period)==1
        sel_mom  = triu(true(size(all_fitted,2)),1);
    else
        % cell to select moments
        sel_mom = cell(length(u_period),1);
        sel_mom(:,1) = {triu(true(size(all_fitted,2)),1)};
    end
      
    % Number of moments and generate matrix to store plug-in estimates
    n_moments = sum(sum(triu(ones(size(dimensions,2),size(dimensions,2)))));

    
    %%% Plugin estimates
    if length(u_period)==1 || isempty(period) 
        cov_mat = cov(all_fitted);   
        plugin = [diag(cov_mat)' cov_mat(sel_mom)'];
    else
        cell_cov_mat = splitapply(@cov2cell,all_fitted,period);
        cell_var = cellfun(@diag,cell_cov_mat,'UniformOutput',false);
        cell_cov = cellfun(@mymoment,cell_cov_mat,sel_mom,'UniformOutput',false);
                
        % Plugin estimates into matrix
        plugin = cell2mat([cellfun(@transpose,cell_var,'UniformOutput',false),...
                                        cellfun(@transpose,cell_cov,'UniformOutput',false)]);             
    end
    
    % Determine label of second order moments
    if size(plugin,2)==3
        str_disp = [{'var_worker'}; {'var_firm'}; {'cov_worker_firm'}];
    else
        tmp_var = strtrim(sprintfc('var_%d ', 3:size(dimensions,2)));
        comb = nchoosek(1:size(dimensions,2),2); % get all the combinations 
        tmp_cov = strtrim(sprintf('cov_%d_%d ', reshape([comb(:,1)'; comb(:,2)'],[],1)));
        str_disp = [{'var_worker'};{'var_firm'}; split(tmp_var)'; split(tmp_cov)];
    end


    disp('Plugin estimates: ')
    T = array2table(plugin);
    T = renamevars(T,1:width(T),str_disp');
    disp(T)
    writetable(T,[filename,'_plugin_estimates.csv'])
   
    %% VARIANCE-COVARIANCE MATRIX ESTIMATION CORRECTIONS
    
    if (strcmp(type_boot,'block') || strcmp(type_boot,'wild_block'))
        group = findgroups(group);
        G = max(group);
        hc = G ./ (G - 1) .* NT ./ dof;
        n_problems_lev = NaN;
        est_var_r = r.^2.*hc;
	elseif (strcmp(type_boot,'diag') && strcmp(type_hc,'hom'))
		est_var_r = NT ./ dof .* (mean(r.^2)).*ones(NT,1);
        n_problems_lev = NaN;
    elseif (strcmp(type_boot,'diag') && strcmp(type_hc,'hc_0'))
        est_var_r = r.^2;
        n_problems_lev = NaN;
    elseif (strcmp(type_boot,'diag') && strcmp(type_hc,'hc_1'))
        hc = NT ./ dof;
        est_var_r = r.^2.*hc;
        n_problems_lev = NaN;
    elseif (strcmp(type_boot,'diag') && strcmp(type_hc,'hc_2') )
        
        %% Estimation of leverage
        fprintf('\n')
        disp('---------- Leverage estimation ---------- ')
        
        % Output is the estimated leverage (hii) and the estimated 1 -
        % leverage (mii)
        tic
        disp('Start to estimate the leverages...')
        [hii, mii, correction_lev] = lev_diagonal_par(n_lev, NT, X, xx, L, mytol,lap);
        disp('Done!')
        toc
        
        %% Diagnostic of leverage estimation
        fprintf('\n')
        disp('---------- Diagnostic ---------- ')
        tic
        [hc_2,n_problems_lev] = diagnostic_par(correction_lev, NT, X, xx, L, mii, mytol,lap);
        toc
        
		hc = hc_2; % without residual to the square
		est_var_r = r.^2.*hc; % this is hc_2 in the paper      
        
        % Check if remaining issues
        if (sum(est_var_r<0)>0)
           disp('THERE ARE REMAINING ISSUES') 
        end
    else 
        disp('Unrecognized type of bootstrap')
    end
    
    %%   BOOTSTRAP CORRECTION 
    fprintf('\n')
    disp('---------- Bootstrap ---------- ')
     
    % Run bootstrap
    tic
    [delta, corrected] = boot_reg_periods_par(n_boot, est_var_r,type_boot,group, design_mat, dimensions, period, sel_mom, n_moments, X, xx, L, plugin, mytol,lap);
    toc
    
    disp('Corrected estimates: ')
    T = array2table(corrected);
    T = renamevars(T,1:width(T),str_disp');
    disp(T)
    writetable(T,[filename,'_corrected_estimates.csv'])

    %% Variance decomposition
    fprintf('\n')
    disp('---------- Variance decomposition ---------- ')
    
    if u_period==1
        var_y = var(y);
    else 
        var_y = splitapply(@var,y,period);
    end
    
    
    % Plugin
    decomp_levels_pi = [var_y plugin(:,1:(size(dimensions,2))) 2.*plugin(:,(size(dimensions,2)+1):end)]; %variance dependent, variance covariates, 2*covariance of covariates, residual variance
    if u_period==1
        decomp_levels_pi = [decomp_levels_pi var_y-sum(decomp_levels_pi(:,2:end))]; % residual variance  
    else
        decomp_levels_pi = [decomp_levels_pi var_y-(splitapply(@(mat) sum(mat),decomp_levels_pi,u_period)-var_y)]; % residual variance per period 
    end
    
    decomp_perc_pi = decomp_levels_pi./decomp_levels_pi(:,1);
    decomp_pi = [u_period decomp_levels_pi; u_period  decomp_perc_pi];
    
    % Determine label of second order moments for the decomposition
    if size(plugin,2)==3
        str_disp = [{'var_worker'}; {'var_firm'}; {'2*cov_worker_firm'}];
    else
        tmp_var = strtrim(sprintfc('var_%d ', 3:size(dimensions,2)));
        comb = nchoosek(1:size(dimensions,2),2); % get all the combinations 
        tmp_cov = strtrim(sprintf('cov_%d_%d ', reshape([comb(:,1)'; comb(:,2)'],[],1)));
        str_disp = [{'var_worker'};{'var_firm'}; split(tmp_var)'; strcat('2*',split(tmp_cov))];
    end

    disp('Variance decomposition of plugin: ')
    str_decomp = [{'period'};{'var_y'};str_disp;{'var_resid'}];
    T = array2table(decomp_pi);
    T = renamevars(T,1:width(T),str_decomp');
    % Add type of decomposition
    type_decomp = split(strtrim([repmat('levels ',1,size(var_y,1)), repmat('percent ',1,size(var_y,1))])); % split converts to cell
    T = addvars(T,type_decomp,'Before',"period");
    disp(T)
    % Store
    writetable(T,[filename,'_var_decomp_plugin.csv'])


    % Corrected
    decomp_levels_b = [var_y corrected(:,1:(size(dimensions,2))) 2.*corrected(:,(size(dimensions,2)+1):end)]; %variance dependent, variance covariates, 2*covariance of covariates, residual variance
    
    if u_period==1
        decomp_levels_b = [decomp_levels_b var_y-sum(decomp_levels_b(:,2:end))]; % residual variance  
    else
        decomp_levels_b = [decomp_levels_b var_y-(splitapply(@(mat) sum(mat),decomp_levels_b,u_period)-var_y)]; % residual variance per period 
    end
    decomp_perc_b = decomp_levels_b./decomp_levels_b(:,1);
    decomp_b = [u_period  decomp_levels_b; u_period decomp_perc_b];
    
    disp('Variance decomposition of corrected: ')
    T = array2table(decomp_b);
    T = renamevars(T,1:width(T),str_decomp');
    % Add type of decomposition
    type_decomp = split(strtrim([repmat('levels ',1,size(var_y,1)), repmat('percent ',1,size(var_y,1))])); % split converts to cell
    T = addvars(T,type_decomp,'Before',"period");
    disp(T)
    writetable(T,[filename,'_var_decomp_corrected.csv'])
    
    

    disp('FINISHED!')
end

