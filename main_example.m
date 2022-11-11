    % ----------------------------------------------------------------------
    %   HOUSEKEEPING
    % ----------------------------------------------------------------------
    clc;
    clear;
    close all;
    warning('off', 'all') 
    diary mydiary_example

	%Set seed
    rng(14);

% % run this once 
%     %%% Install 
%     %Install CMG Routine: downloaded from here: http://www.cs.cmu.edu/~jkoutis/cmg.html  
%     cd CMG;    
%     path(path,'CMG'); %this contains the main LeaveOut Routines.
%     MakeCMG;
%     cd ..
% %

    path(path,'matlab_bgl/'); %note: matlab BGL obtained from http://www.mathworks.com/matlabcentral/fileexchange/10922
    path(path,'src'); %this contains the functions for bootstrap correction.

    % ----------------------------------------------------------------------
    %   SOURCES AND PARAMETERS
    % ----------------------------------------------------------------------   
    
    % Bootstrap parameters    
    N_boot = 300; 
    type_boot = 'diag';
    type_hc = 'hc_2';
    group = [];
    period = [];
    n_boot = N_boot;
    n_lev = N_boot;
    mytol = 1e-6;
    
      
    % Source data
    nm_data = '.\data\example.csv'; 
    


    % ----------------------------------------------------------------------
    %   LOAD DATA
    % ----------------------------------------------------------------------
    
    %% READ DATA
    data=readtable(nm_data);
    
    id=data.id;
    firmid=data.firmid;
    occ_id=data.occ_id;
    pe_t_0=data.pe_t;
    fe_t_0=data.fe_t;
    occ_t=data.occ_t;
    age = data.age;
    age_sq = age.^2;
    y = data.y;
    year = data.year;

   
    clear data


    % Fixed effects, controls
    mat_controls = [age age_sq];
    mat_id_fe = [id, firmid, occ_id];

    % ----------------------------------------------------------------------
    %   A)
    % ----------------------------------------------------------------------    
    disp('------------------------------------')
    disp('  EXAMPLE RESIDUALIZING')
    disp('------------------------------------')
    
    try
        pool=parpool;
    end


    % Additional parameters
    resid_controls = 1;
    nm_out = '.\results\example_resid';
	

    start=tic;
    [plugin,delta,corrected,decomp_pi,decomp_b,dimensions,NT,n_problems_lev] = boot_correction(y,mat_id_fe,mat_controls,n_lev,n_boot,period,group,type_boot,type_hc,mytol,resid_controls,nm_out);
    time= toc(start);       

    
    %% Store results    
    save([nm_out,'.mat'],'plugin','corrected','decomp_pi','decomp_b','time','type_boot','type_hc','n_boot','n_lev','n_problems_lev')


    delete(gcp)

    % ----------------------------------------------------------------------
    %   B)
    % ----------------------------------------------------------------------

    disp('------------------------------------')
    disp('  EXAMPLE WITHOUT RESIDUALIZING')
    disp('------------------------------------')
    
    try
        pool=parpool;
    end

    % Fixed effects, controls
    mat_controls = [age age_sq];
    mat_id_fe = [id, firmid, occ_id];

    
    %% PARAMETERS
   % Additional parameters
    resid_controls = 0;
    nm_out = '.\results\example_no_resid';

    
    %% Bootstrap corrections

    start=tic;
    [plugin,delta,corrected,decomp_pi,decomp_b,dimensions,NT,n_problems_lev] = boot_correction(y,mat_id_fe,mat_controls,n_lev,n_boot,period,group,type_boot,type_hc,mytol,resid_controls,nm_out);
    time= toc(start);

  
        
    %% Store results    
    save([nm_out,'.mat'],'plugin','corrected','decomp_pi','decomp_b','time','type_boot','type_hc','n_boot','n_lev','nm_data','n_problems_lev')

    delete(gcp)


    % ----------------------------------------------------------------------
    %   C)
    % ----------------------------------------------------------------------    
    disp('------------------------------------')
    disp('  EXAMPLE RESIDUALIZING. Different periods')
    disp('------------------------------------')
    
    try
        pool=parpool;
    end


    % Additional parameters
    resid_controls = 1;
    nm_out = '.\results\example_resid_py';
    period = year;
	

    start=tic;
    [plugin,delta,corrected,decomp_pi,decomp_b,dimensions,NT,n_problems_lev] = boot_correction(y,mat_id_fe,mat_controls,n_lev,n_boot,period,group,type_boot,type_hc,mytol,resid_controls,nm_out);
    time= toc(start);       

    
    %% Store results    
    save([nm_out,'.mat'],'plugin','corrected','decomp_pi','decomp_b','time','type_boot','type_hc','n_boot','n_lev','n_problems_lev')

    delete(gcp)

    diary off
 
   