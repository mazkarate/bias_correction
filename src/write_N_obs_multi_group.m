function write_N_obs_multi_group(ind_export, u_group_old_sorted,...
    N_multi_group, N_obs_group, N_id_group, N_firmid_group, filename, v_filename_group)

    for i= 1:N_multi_group
 
            if ind_export == 1 && size(v_filename_group,2) ==0
                writetable(array2table([["overall"; u_group_old_sorted{i}] N_obs_group{i} N_id_group{i} N_firmid_group{i}]), strcat(filename,'_g',num2str(i),'_N_obs_group.csv'))
            elseif ind_export == 1
                writetable(array2table([["overall"; u_group_old_sorted{i}] N_obs_group{i} N_id_group{i} N_firmid_group{i}]), strcat(filename,'_',v_filename_group{i},'_N_obs_group.csv'))
            end
    
    end

end
