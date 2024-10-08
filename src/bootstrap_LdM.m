function [delta, corrected] = bootstrap_LdM(n_boot, B_plus, B_minus, dim_fe_mover, X, xx, L,...
                                        plugin, mytol, NT_mover, id, firmid, firmid_stayers,...
                                        N_obs_per_firmid_y, firmid_y, F_y)

    % Pre-assign estimate of bias
    delta = 0;
    
    NT = size(B_plus,1); % Change length of simulating vectors
    
    % B is a vector
    if ~isempty(B_minus)
        parfor s=1:n_boot
    
                % Determine the dependent variable for positive and negative
                % eigenvalues
                v_tilde = randsample([-1;1],NT,true);
    
                % Solve normal equations with preconditioning matrix L
                % B+
                y_b = B_plus.*v_tilde;
                [b_tilde,~]    =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
                [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
                
                %%% Average worker fixed effects at the firm per year
                sum_pe = F_y * pe_tilde;    
                sum_pe = sum_pe(firmid_y); % repeat vector
                
                % Average of estimated worker fixed effects for coworkers
                av_coworkers = (sum_pe - pe_tilde)./(N_obs_per_firmid_y-1); % assign average coworker FE by removing the current worker fixed effect for firms with more than 2 workers
                
                aux = compute_moments_LdM(pe_tilde, fe_tilde, av_coworkers);
                
                % B-
                y_b = B_minus.*v_tilde;
                [b_tilde,~]   =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
                [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
                
                %%% Average worker fixed effects at the firm per year
                sum_pe = F_y * pe_tilde;    
                sum_pe = sum_pe(firmid_y); % repeat vector
                
                % Average of estimated worker fixed effects for coworkers
                av_coworkers = (sum_pe - pe_tilde)./(N_obs_per_firmid_y-1); % assign average coworker FE by removing the current worker fixed effect for firms with more than 2 workers

                % Delta                        
                delta = delta + (aux - compute_moments_LdM(pe_tilde, fe_tilde, av_coworkers));
    
        end
    else
        parfor s=1:n_boot
    
                % Determine the dependent variable for only positive eigenvalues
                v_tilde = randsample([-1;1],NT,true);
    
                % Solve normal equations with preconditioning matrix L
                % B+
                y_b = B_plus.*v_tilde;
                [b_tilde,~]    =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
                [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
                
                %%% Average worker fixed effects at the firm per year
                sum_pe = F_y * pe_tilde;    
                sum_pe = sum_pe(firmid_y); % repeat vector
                
                % Average of estimated worker fixed effects for coworkers
                av_coworkers = (sum_pe - pe_tilde)./(N_obs_per_firmid_y-1); % assign average coworker FE by removing the current worker fixed effect for firms with more than 2 workers

                    
                % Delta                        
                delta = delta + (compute_moments_LdM(pe_tilde, fe_tilde, av_coworkers));
    
        end
    end  

    %Store bootstrap results
    delta = delta./n_boot; % take the average
	corrected = plugin - delta;

end