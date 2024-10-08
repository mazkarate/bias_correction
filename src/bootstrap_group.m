function [delta, corrected] = bootstrap_group(n_boot, B_plus, B_minus, dim_fe_mover, X, xx, L, plugin, mytol, NT_mover,...
                                                                       id, firmid, firmid_stayers, i_group, ch_group)
                                 
    % Pre-assign estimate of bias
    delta = 0;
    
    NT = size(B_plus,1); %Change length of simulating vectors
    
    % B is a vector
    if ~isempty(B_minus)
        parfor s=1:n_boot
    
                % Determine the dependent variable for positive and negative
                % eigenvalues
                v_tilde = randsample([-1;1],NT,true);
    
                % run AKM with preconditioning matrix L
                % B+
                y_b = B_plus.*v_tilde;
                [b_tilde,~]    =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);                
                [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
                aux = compute_moments_group(pe_tilde, fe_tilde, i_group, ch_group);
                
                % B-
                y_b = B_minus.*v_tilde;
                [b_tilde,~]   =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
                [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
                aux = aux - compute_moments_group(pe_tilde, fe_tilde, i_group, ch_group);
    
                % Delta                        
                delta = delta + aux;
    
        end
    else
        parfor s=1:n_boot
    
                % Determine the dependent variable for only positive eigenvalues
                v_tilde = randsample([-1;1],NT,true);
    
                % run AKM with preconditioning matrix L
                % B+
                y_b = B_plus.*v_tilde;
                [b_tilde,~]    =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
                [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
                aux = compute_moments_group(pe_tilde, fe_tilde, i_group, ch_group);
                    
                % Delta                        
                delta = delta + aux;
    
        end
    end  

    %Store bootstrap results
    delta = delta./n_boot; % take the average
	corrected = plugin - delta;

end