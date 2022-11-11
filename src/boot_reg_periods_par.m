function [delta, corrected] = boot_reg_periods_par(n_boot, est_var_r,type_boot,group, design_mat, dimensions, period, sel_mom, n_moments, X, xx, L, plugin, mytol, lap)
    

 
    NT = length(est_var_r);
   
    if exist('period')==0 || isempty(period)
       u_period = 1;
    else
       u_period = unique(period);
    end
    
    % Pre-assign estimate of bias
    if length(u_period)==1 || isempty(period)
        delta_s = zeros(n_boot,n_moments);  
    else
        delta_s = cell(length(u_period),1); % each element in the cell is a bootstrap
        % Initialize cell with empty matrices
        delta_s = {zeros(n_boot,n_moments)};
    end

    % Pass variables to workers
    X_par = parallel.pool.Constant(X);
    xx_par = parallel.pool.Constant(xx);
    design_mat_par = parallel.pool.Constant(design_mat);
    est_var_r_par = parallel.pool.Constant(est_var_r);

    if lap ==1

    % Pass preconditioner to workers
    L_par = parallel.pool.Constant(@() L);

        % Determine cases
        if (strcmp(type_boot,'diag'))
            
            
            parfor s=1:n_boot
                    
                    % Determine the dependent variable
                    v_tilde = randsample([-1;1],NT,true);
                    y_b = sqrt(est_var_r_par.Value).*v_tilde;
                                       
                    % run AKM with preconditioning matrix L and initial guess b
                    [b_tilde,~]=pcg(xx_par.Value,X_par.Value'*y_b,mytol,300,L_par.Value);
                    
                    % cell with betas
                    betas = cell(1,size(dimensions,2));
                    for i=1:size(dimensions,2)
                        if i==1
                            pos = 0;
                        else
                            pos = pos + dimensions(i-1);
                        end
                        betas{i} = b_tilde((pos+1):(pos+dimensions(i)));
                    end
                    
                    cell_all_fitted = cellfun(@(x,y) x*y, design_mat_par.Value,betas, 'UniformOutput',false);
                    
                    % One fitted per column
                    all_fitted = reshape(cell2mat(cell_all_fitted),NT,[]);
                    
                    %%% Delta
                    if length(u_period)==1 || isempty(period)
                        cov_mat = cov(all_fitted);   
                        delta_s(s,:) = [diag(cov_mat)' cov_mat(sel_mom)'];
                    else                 
						cell_cov_mat = splitapply(@cov2cell,all_fitted,period);
						cell_var = cellfun(@diag,cell_cov_mat,'UniformOutput',false);
						cell_cov = cellfun(@mymoment,cell_cov_mat,sel_mom,'UniformOutput',false);
							
						% Arrange to get variances and covariances in one line
						% vector
						delta_s{s,:} = cell2mat([cellfun(@transpose,cell_var,'UniformOutput',false),...
												 cellfun(@transpose,cell_cov,'UniformOutput',false)]);
                    end
                    
            end
            
        elseif (strcmp(type_boot,'wild_block')) 
                G = max(group);
             parfor s=1:n_boot
                    
                    % Determine the dependent variable
                    v_tilde = randsample([-1;1],G,true); % Random variable per group
                    y_b = sqrt(est_var_r).*v_tilde(group); % evaluate the Rademacher at the group level
                    
                    % run AKM with preconditioning matrix L and initial guess b
                    [b_tilde,~]=pcg(xx_par.Value,X_par.Value'*y_b,mytol,300,L_par.Value);
                    
                    % cell with betas
                    betas = cell(1,size(dimensions,2));
                    for i=1:size(dimensions,2)
                        if i==1
                            pos = 0;
                        else
                            pos = pos + dimensions(i-1);
                        end
                        betas{i} = b_tilde((pos+1):(pos+dimensions(i)));
                    end
                    
                    cell_all_fitted = cellfun(@(x,y) x*y, design_mat_par.Value,betas, 'UniformOutput',false);
                    
                    % One fitted per column
                    all_fitted = reshape(cell2mat(cell_all_fitted),NT,[]);
                    
                    %%% Delta
                    if length(u_period)==1 || isempty(period)
                        cov_mat = cov(all_fitted);   
                        delta_s(s,:) = [diag(cov_mat)' cov_mat(sel_mom)'];
                    else    
                        cell_cov_mat = splitapply(@cov2cell,all_fitted,period);
                        cell_var = cellfun(@diag,cell_cov_mat,'UniformOutput',false);
                        cell_cov = cellfun(@mymoment,cell_cov_mat,sel_mom,'UniformOutput',false);
                        
                       % Arrange to get variances and covariances in one line
                       % vector
                        delta_s{s,:} = cell2mat([cellfun(@transpose,cell_var,'UniformOutput',false),...
                                                cellfun(@transpose,cell_cov,'UniformOutput',false)]);
                    end
                    
            end
            
        elseif (strcmp(type_boot,'block'))   
                        disp('Block bootstrap not supported yet.')
        
                        return
                        
        end  
            %   Store bootstrap results
            if length(u_period)==1 || isempty(period)
                delta = mean(delta_s,1); % mean by columns
                corrected = plugin - delta; 
            else 
                delta = mean(cat(3,delta_s{:}),3); % take the average across cell bins
                corrected = plugin - delta;
            end
    
    else
       if size(L,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
            % Pass preconditioner to workers
            L_par = parallel.pool.Constant(L);
       end
       
        % Determine cases
        if (strcmp(type_boot,'diag'))
            
            
            parfor s=1:n_boot
                    
                    % Determine the dependent variable
                    v_tilde = randsample([-1;1],NT,true);
                    y_b = sqrt(est_var_r_par.Value).*v_tilde;
                            
                    % run AKM with preconditioning matrix L and initial guess b
                    if size(L,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
                            [b_tilde,~]=pcg(xx_par.Value,X_par.Value'*y_b,mytol,300,L_par.Value,L_par.Value');
                    else 
                            [b_tilde,~]=pcg(xx_par.Value,X_par.Value'*y_b,mytol,300);
                    end
                    
                    % cell with betas
                    betas = cell(1,size(dimensions,2));
                    for i=1:size(dimensions,2)
                        if i==1
                            pos = 0;
                        else
                            pos = pos + dimensions(i-1);
                        end
                        betas{i} = b_tilde((pos+1):(pos+dimensions(i)));
                    end
                    
                    cell_all_fitted = cellfun(@(x,y) x*y, design_mat_par.Value,betas, 'UniformOutput',false);
                    
                    % One fitted per column
                    all_fitted = reshape(cell2mat(cell_all_fitted),NT,[]);
                    
                    %%% Delta
                    if length(u_period)==1 || isempty(period)
                        cov_mat = cov(all_fitted);   
                        delta_s(s,:) = [diag(cov_mat)' cov_mat(sel_mom)'];
                    else       
						cell_cov_mat = splitapply(@cov2cell,all_fitted,period);
						cell_var = cellfun(@diag,cell_cov_mat,'UniformOutput',false);
						cell_cov = cellfun(@mymoment,cell_cov_mat,sel_mom,'UniformOutput',false);
							
						% Arrange to get variances and covariances in one line
						% vector
						delta_s{s,:} = cell2mat([cellfun(@transpose,cell_var,'UniformOutput',false),...
												 cellfun(@transpose,cell_cov,'UniformOutput',false)]);
                    end
                    
            end
            
        elseif (strcmp(type_boot,'wild_block')) 
                G = max(group);
             parfor s=1:n_boot
                    
                    % Determine the dependent variable
                    v_tilde = randsample([-1;1],G,true); % Random variable per group
                    y_b = sqrt(est_var_r).*v_tilde(group); % evaluate the Rademacher at the group level
                    
                    % run AKM with preconditioning matrix L and initial guess b
                    if size(L,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
                            [b_tilde,~]=pcg(xx_par.Value,X_par.Value'*y_b,mytol,300,L_par.Value,L_par.Value');
                    else 
                            [b_tilde,~]=pcg(xx_par.Value,X_par.Value'*y_b,mytol,300);
                    end
                    
                    % cell with betas
                    betas = cell(1,size(dimensions,2));
                    for i=1:size(dimensions,2)
                        if i==1
                            pos = 0;
                        else
                            pos = pos + dimensions(i-1);
                        end
                        betas{i} = b_tilde((pos+1):(pos+dimensions(i)));
                    end
                    
                    cell_all_fitted = cellfun(@(x,y) x*y, design_mat_par.Value,betas, 'UniformOutput',false);
                    
                    % One fitted per column
                    all_fitted = reshape(cell2mat(cell_all_fitted),NT,[]);
                    
                    %%% Delta
                    if length(u_period)==1 || isempty(period)
                        cov_mat = cov(all_fitted);   
                        delta_s(s,:) = [diag(cov_mat)' cov_mat(sel_mom)'];
                    else
                        
                        cell_cov_mat = splitapply(@cov2cell,all_fitted,period);
                        cell_var = cellfun(@diag,cell_cov_mat,'UniformOutput',false);
                        cell_cov = cellfun(@mymoment,cell_cov_mat,sel_mom,'UniformOutput',false);
                        
                       % Arrange to get variances and covariances in one line
                       % vector
                        delta_s{s,:} = cell2mat([cellfun(@transpose,cell_var,'UniformOutput',false),...
                                                cellfun(@transpose,cell_cov,'UniformOutput',false)]);
                    end
                    
            end
            
        elseif (strcmp(type_boot,'block'))
                        disp('Block bootstrap not supported yet.')
        
                        return
                        
        end
        
		%   Store bootstrap results
		if length(u_period)==1 || isempty(period)
			delta = mean(delta_s,1); % mean by columns
			corrected = plugin - delta; 
		else 
			delta = mean(cat(3,delta_s{:}),3); % take the average across cell bins
			corrected = plugin - delta;
		end

    end
end