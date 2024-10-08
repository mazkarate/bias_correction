function [est_var_r_collapsed, weight_stayers] = collapse_var_stayers(est_var_r, dim_fe_mover, id)

    % Indicator of stayer (last ids) 
    ind_stayer = id > dim_fe_mover(1);
    % Id's of stayers starting from 1
    id_stayer_norm = id(ind_stayer) - dim_fe_mover(1);
        
    % Get the number of periods that stayers are in the sample
    weight_stayers = accumarray(id_stayer_norm, 1); % Find number of counts per stayer in sample

    est_var_r_collapsed_stayers = accumarray(id_stayer_norm, est_var_r(ind_stayer), [], @(x) sum(x)); % Sum variances of obs per stayer
    est_var_r_collapsed_stayers = est_var_r_collapsed_stayers./(weight_stayers.^2); % Divide by time in sample, squared

    % Export the estimated variances for movers (one per observations) and
    % for stayers (one per match)
    est_var_r_collapsed = [est_var_r(~ind_stayer); est_var_r_collapsed_stayers];