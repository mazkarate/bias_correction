function [delta, corrected] = bootstrap(n_boot, B_plus, B_minus, dim_fe_mover, X, xx, L,...
                                        plugin, mytol, NT_mover, id, firmid, firmid_stayers)

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
                aux = compute_moments(pe_tilde, fe_tilde);
                
                % B-
                y_b = B_minus.*v_tilde;
                [b_tilde,~]   =   pcg(xx, X'*y_b(1:NT_mover), mytol, 300, L);
                [pe_tilde, fe_tilde] = compute_fe_boot(b_tilde, dim_fe_mover, y_b(NT_mover + 1 : end), id, firmid, firmid_stayers);
                aux = aux - compute_moments(pe_tilde, fe_tilde);
    
                % Delta                        
                delta = delta + aux;
    
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
                aux = compute_moments(pe_tilde, fe_tilde);
                    
                % Delta                        
                delta = delta + aux;
    
        end
    end  

    %Store bootstrap results
    delta = delta./n_boot; % take the average
	corrected = plugin - delta;

end