function write_est_multi_group(ind_export, est, nm_est, str_disp, u_group_old_sorted, N_multi_group, N_obs_group, N_id_group, N_firmid_group, filename, v_filename_group)

    for i= 1:N_multi_group
    
            disp([nm_est,' estimates: ']) 
    
            disp(['Grouping: ', num2str(i)]) 
            T = array2table([est{i} N_obs_group{i} N_id_group{i} N_firmid_group{i}]);
            T = renamevars(T,1:width(T),str_disp');
            T = addvars(T,["overall"; u_group_old_sorted{i}],'Before',"var_worker");
            T = renamevars(T,"Var1","group_old");          
            % Display only 20 lines
            if size(est{i},1)<20
                disp(T)
            else 
                disp(T(1:20,:))
            end
            if ind_export == 1 && size(v_filename_group,2) ==0
                writetable(T, strcat(filename,'_g',num2str(i),'_',nm_est,'_estimates.csv'))   
            elseif ind_export == 1
                writetable(T, strcat(filename,'_',v_filename_group{i},'_',nm_est,'_estimates.csv'))
            end
    
    end

end
