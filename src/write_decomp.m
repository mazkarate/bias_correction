function [decomp_pi, decomp_b] = write_decomp(ind_export, plugin, corrected, str_disp, var_y,...
        dim_fe, filename)
       
        % Plugin
        decomp_levels_pi = [var_y plugin(:,1:(size(dim_fe,2)))...
            2.*plugin(:,(size(dim_fe,2)+1):end)];
        decomp_levels_pi = [decomp_levels_pi var_y - ...
            (splitapply(@(mat) sum(mat),decomp_levels_pi,(1:size(decomp_levels_pi,1))')-var_y)]; % Residual variance per group 
        
        decomp_perc_pi = decomp_levels_pi./decomp_levels_pi(:,1);
        decomp_pi = [decomp_levels_pi; decomp_perc_pi];   
    
        disp('Variance decomposition of plugin: ')
        str_decomp = [{'var_y'};str_disp;{'var_resid'}];
        T = array2table(decomp_pi);
        T = renamevars(T,1:width(T),str_decomp');
        % Add type of decomposition
        type_decomp = split(strtrim([repmat('levels ',1,size(var_y,1)), repmat('percent ',1,size(var_y,1))])); % Split converts to cell
        T = addvars(T,type_decomp,'Before',"var_y");
        disp(T)
        
        if ind_export == 1
            writetable(T, strcat(filename,'_var_decomp_plugin.csv'))   
        end

        % Corrected
        decomp_levels_b = [var_y corrected(:,1:(size(dim_fe,2)))...
            2.*corrected(:,(size(dim_fe,2)+1):end)];    
        decomp_levels_b = [decomp_levels_b var_y - ...
            (splitapply(@(mat) sum(mat),decomp_levels_b,(1:size(decomp_levels_b,1))')-var_y)]; % Residual variance per group 
        decomp_perc_b = decomp_levels_b./decomp_levels_b(:,1);
        decomp_b = [decomp_levels_b; decomp_perc_b];   
    
        disp('Variance decomposition of corrected: ')
        T = array2table(decomp_b);
        T = renamevars(T,1:width(T),str_decomp');
        % Add type of decomposition
        type_decomp = split(strtrim([repmat('levels ',1,size(var_y,1)), repmat('percent ',1,size(var_y,1))])); % Split converts to cell
        T = addvars(T,type_decomp,'Before',"var_y");
        disp(T)
        
        if ind_export == 1
            writetable(T, strcat(filename,'_var_decomp_corrected.csv'))   
        end
end
