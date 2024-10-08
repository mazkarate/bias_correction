function [decomp_pi, decomp_b] = write_decomp_multi_group(ind_export, plugin, corrected, str_disp, var_y,...
    u_group_old_sorted, dim_fe, N_multi_group, filename, v_filename_group)

    for i=1:N_multi_group

        disp(['Grouping: ', num2str(i)]) 
        
        % Plugin
        decomp_levels_pi = [var_y{i} plugin{i}(:,1:(size(dim_fe,2)))...
            2.*plugin{i}(:,(size(dim_fe,2)+1):end)];
        decomp_levels_pi = [decomp_levels_pi var_y{i} - ...
            (splitapply(@(mat) sum(mat),decomp_levels_pi,(1:size(decomp_levels_pi,1))')-var_y{i})]; % Residual variance per group 
        
        decomp_perc_pi = decomp_levels_pi./decomp_levels_pi(:,1);
        decomp_pi = [decomp_levels_pi; decomp_perc_pi];   
    
        disp('Variance decomposition of plugin: ')
        str_decomp = [{'var_y'};str_disp;{'var_resid'}];
        T = array2table(decomp_pi);
        T = renamevars(T,1:width(T),str_decomp');
        % Add type of decomposition
        type_decomp = split(strtrim([repmat('levels ',1,size(var_y{i},1)), repmat('percent ',1,size(var_y{i},1))])); % Split converts to cell
        T = addvars(T,type_decomp,["overall";u_group_old_sorted{i} ; "overall"; u_group_old_sorted{i} ],'Before',"var_y");
        T = renamevars(T,"Var2","group_old");

        if size(plugin{i},1)<20
            disp(T)
        else 
            disp(T(1:20,:))
        end
        if ind_export == 1 && size(v_filename_group,2) ==0
            writetable(T, strcat(filename,'_g',num2str(i),'_var_decomp_plugin.csv'))   
        elseif ind_export == 1
            writetable(T, strcat(filename,'_',v_filename_group{i},'_var_decomp_plugin.csv'))
        end

        % Corrected
        decomp_levels_b = [var_y{i} corrected{i}(:,1:(size(dim_fe,2)))...
            2.*corrected{i}(:,(size(dim_fe,2)+1):end)];    
        decomp_levels_b = [decomp_levels_b var_y{i} - ...
            (splitapply(@(mat) sum(mat),decomp_levels_b,(1:size(decomp_levels_b,1))')-var_y{i})]; % Residual variance per group 
        decomp_perc_b = decomp_levels_b./decomp_levels_b(:,1);
        decomp_b = [decomp_levels_b; decomp_perc_b];   
    
        disp('Variance decomposition of corrected: ')
        T = array2table(decomp_b);
        T = renamevars(T,1:width(T),str_decomp');
        % Add type of decomposition
        T = addvars(T,type_decomp,["overall";u_group_old_sorted{i} ; "overall";u_group_old_sorted{i} ],'Before',"var_y");
        T = renamevars(T,"Var2","group_old");
        % Display only 20 lines
        if size(plugin{i},1)<20
            disp(T)
        else 
            disp(T(1:20,:))
        end
        if ind_export == 1 && size(v_filename_group,2) ==0
            writetable(T, strcat(filename,'_g',num2str(i),'_var_decomp_corrected.csv'))   
        elseif ind_export == 1
            writetable(T, strcat(filename,'_',v_filename_group{i},'_var_decomp_corrected.csv'))
        end
    end

end
