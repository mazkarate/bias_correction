function write_est(ind_export, est, nm_est, str_disp, N_obs_group, N_id_group, N_firmid_group, filename)
    
    disp([nm_est,' estimates: ']) 
    T = array2table([est N_obs_group N_id_group N_firmid_group]);
    T = renamevars(T,1:width(T),str_disp');
    disp(T)

    if ind_export == 1
        writetable(T, strcat(filename,'_',nm_est,'_estimates.csv'))   
    end



end