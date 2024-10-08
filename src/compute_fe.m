function [pe_tilde, fe_tilde] = compute_fe(b, dim_fe_mover, y, id, firmid, NT_mover)

    % Indicator of stayer (last ids) 
    ind_stayer = id > dim_fe_mover(1);
    % Id's of stayers starting from 1
    id_stayer_norm = id(ind_stayer) - dim_fe_mover(1); 
    
    % Value firm FE
    b_firm = [-b( ( dim_fe_mover(1) + 1 ):end ); 0]; % Need to add zero as last firm fixed effect is normalized. Added "minus" because of Laplacian representation 

    % Estimated FE for movers
    pe_tilde_mover = b( id( ~ind_stayer ) );
    fe_tilde_mover = b_firm( firmid( ~ind_stayer ) );

    % Estimated FE for stayers
    fe_tilde_stayer = b_firm( firmid( ind_stayer ) ); % Stayers in the connected set must be at a firm of a mover

    % Firm fixed effects
    fe_tilde = [fe_tilde_mover; fe_tilde_stayer];

    % Person fixed effects for stayers
    pe_tilde_stayer = accumarray(id_stayer_norm, y( NT_mover + 1:end ), [], @(x) mean(x))...
                        - accumarray(id_stayer_norm, fe_tilde_stayer ,[], @(x) mean(x)); % One estimate per id

    pe_tilde_stayer = pe_tilde_stayer( id_stayer_norm ); % Repeat estimates for number of observations of each stayer

    % All estimated person FE (movers, then stayers)   
    pe_tilde = [pe_tilde_mover; pe_tilde_stayer];

end