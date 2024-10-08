function [pe_tilde, fe_tilde] = compute_fe_boot(b, dimensions_mover, y, id, firmid, firmid_stayers)

    %Number of movers (number of worker fixed effects for movers)
    N_movers = dimensions_mover(1);
    
    %% Get vector of Firm Fixed Effcts

    %%% Value firm FE estimated
    b_firm = [-b( N_movers + 1 : end ); 0]; %Need to add zero as last firm fixed effect is normalized. Added "minus" because of Laplacian representation
    
    %%% Vector of Firm fixed effects
    fe_tilde = b_firm(firmid);
    
    %% Get vector of Worker Fixed Effects

    % Estimated Firm FE for stayers by id
    fe_tilde_stayer = b_firm(firmid_stayers); % stayers in the connected set must be at a firm of a mover
    
    % Estimated person fixed effects for stayers by id
    b_stayer = y - fe_tilde_stayer;
    
    % Estimated worker fixed effects for both movers and stayers
    b_worker = [b( 1 : N_movers ); b_stayer];
    
    %%% Vector of Worker fixed effects (movers then stayers)   
    pe_tilde = b_worker(id);
end