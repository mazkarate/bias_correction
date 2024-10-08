function [true_moments, NT] = compute_true_moments_cluster(id, firmid, pe_t_0, fe_t_0, leave_out_level, cluster, year)
    %Code to compute true moments in data. Useful when doing simulations
    
    disp(['Size initial data: ',num2str(length(id))])
    
    %Lagfirmid
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker

    % Rename FE identifiers
    disp('Relabeling ids...')
    NT = length(id);
    sel=~isnan(lagfirmid);

    %relabel the firms
    [~,~,n]=unique([firmid;lagfirmid(sel)]);

    firmid=n(1:NT);
    lagfirmid(sel)=n(NT+1:end);
    
    %relabel the workers
    [~,~,n]=unique(id);
    id=n;
    
    %% Compute the connected set and relabel
    disp('---------- Connected set ---------- ')
    %Find the connected set. Then define dimensions and relevant matrices.
    sel = connected(firmid,lagfirmid,sel);

    %%% Filter after finding connected set
    firmid=firmid(sel); id=id(sel); 
    lagfirmid=lagfirmid(sel);
    if ~isempty(cluster)
        cluster = cluster(sel);
    end
    if ~isempty(year)
        year = year(sel);
    end

    fe_t_0_aux = fe_t_0(sel); pe_t_0_aux = pe_t_0(sel);
    
    % Rename FE ids
    disp('Relabeling ids again...')
    NT = size(id,1);
    sel=~isnan(lagfirmid);
    
    %relabel the firms
    [~,~,n]=unique([firmid; lagfirmid(sel)]);

    firmid=n(1:NT);
    lagfirmid(sel)=n(NT+1:end);

    %relabel the workers
    [~,~,n]=unique(id);
    id=n;
        
    disp(['Size connected set data: ',num2str(size(fe_t_0_aux,1))])
    
    if (strcmp(leave_out_level, 'obs') || strcmp(leave_out_level, 'match') || strcmp(leave_out_level, 'worker'))
        disp('---------- Pruning ---------- ')   
        tic
            % Remove ids that appear only once. They dont generate connections
            % and have leverage equal to 1

            % Count obs per id
            [~,~,iunique] = unique(id);
            n = accumarray(iunique(:),1);

            % Remove observations with only one observation
            sel = (n>1); 
            sel = sel(iunique); % Repeat to filter original vector
            clear iunique
            
            % Filter out one-time ids
            firmid=firmid(sel); id=id(sel); 
            
            disp(['Size after removing one-timers: ',num2str(size(firmid,1))])
              
            % Rename FE ids after removing one-timers
            sel=~isnan(lagfirmid);
    
            % Relabel the firms
            [~,~,n]=unique([firmid;lagfirmid(sel)]);
            firmid=n(1:NT);
                
            % Relabel the workers
            [~,~,id]=unique(id,'stable');

            % Choose pruning type: leave worker, match or observation out
            % Leave observation-out is the default
            if isempty(year)
                disp('No LdM')
                if (strcmp(leave_out_level,'worker')) 
                        sel = pruning_worker(id, firmid); % leave worker as KSS
                elseif (strcmp(leave_out_level,'match')) 
                        sel = pruning_match(id, firmid); % leave match out
                else  
                        sel = pruning_obs(id, firmid); % leave observation out
                end
            else 
                disp('With LdM')
                if (strcmp(leave_out_level,'worker')) 
                        sel = pruning_worker_LdM(id, firmid, year); % leave worker as KSS
                elseif (strcmp(leave_out_level,'match')) 
                        sel = pruning_match_LdM(id, firmid, year); % leave match out
                else  
                        sel = pruning_obs_LdM(id, firmid, year); % leave observation out
                end

            end

            % Filter variables and rename identifiers
            fe_t_0_aux = fe_t_0_aux(sel); pe_t_0_aux = pe_t_0_aux(sel); 

            disp(['Size pruned data: ',num2str(size(fe_t_0_aux,1))])
        toc
        
    elseif strcmp(leave_out_level, 'cluster')
    
    
        fprintf('\n')
        disp('---------- Pruning ---------- ')
        tic
    
        % Check what type of pruning to do: using tripartite graph or looping
        % over the bipartite graph
        firm_worker_cluster = check_type_cluster(id, firmid, cluster);
    
        % Do pruning at the cluster level: selects surviving firms and clusters
        if isempty(year)
            disp('No LdM')
            if firm_worker_cluster
                sel = pruning_cluster(id, firmid, cluster);
            else
                % New cluster vector gets cluster at observation level for those
                % belonging to bad clusters
                [sel, ~] = pruning_cluster_loop(id, firmid, cluster);        
                if isempty(sel)
                    fprintf('\n')
                    error('ERROR: Sample is not leave-cluster-out estimable. Consider different level of clustering.' );            
                end
            end

        else 
            disp('With LdM')
            if firm_worker_cluster
                sel = pruning_cluster_LdM(id, firmid, cluster, year);
            else
                % New cluster vector gets cluster at observation level for those
                % belonging to bad clusters
                [sel, ~] = pruning_cluster_loop_LdM(id, firmid, cluster, year);        
                if isempty(sel)
                    fprintf('\n')
                    error('ERROR: Sample is not leave-cluster-out estimable. Consider different level of clustering.' );            
                end
            end
        end
        % Filter variables
            fe_t_0_aux = fe_t_0_aux(sel); pe_t_0_aux = pe_t_0_aux(sel);      

            disp(['Size pruned data: ',num2str(size(fe_t_0_aux,1))])
        toc
    end
    clear lagfirmid sel
      
    NT = length(fe_t_0_aux);

    true_moments = compute_moments(pe_t_0_aux, fe_t_0_aux);
    
    
