function [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_all(y, id, firmid, other_fe, mat_controls, group,...
			n_lev, n_boot, type_hc, type_leave,...
			cluster, ind_light, mytol, LdM_mom, year,...
			filename, ind_export, v_filename_group)

%{
  
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
    %% 					GENERAL DESCRIPTION
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-  

    This function computes the bootstrap correction of second order
    moments of two-way fixed effects (e.g. in the labor market
    application, a decomposition of log wages into worker and firm fixed
    effects). We use AKM jargon (workers, firms) when describing the code for
    simplicity.
    
    The mandatory input is a person-year dataset that has to be sorted
    by workers' identifiers (id) and year. The function requires as input
    the log wages or the outcome variable and vectors of worker and firm 
    identifiers. The rest of the parameters take default values if they are not provided.

    The function automatically performs the computation of the largest connected set or leave
    out connected set depending on the case. The function computes corrected variance components and a
    variance decomposition.
    
    % Version:
    1.0: First version. 27/09/2020.
    2.0: Replication package. 04/10/2022.
    2.1: Change preconditioner for regressions with two leading fixed effects.
    3.1: 2 step estimation of worker fixed effects for stayers collapsing
    observations. 18/06/2024
    
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
    %% 					DESCRIPTION OF THE INPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                            %-MANDATORY INPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-

    . y: Outcome. Dimensions: N* x 1; N*= # of person-year observations.

    . id: Matrix with worker fixed effect identifiers. Dimensions: N* x 1
     
    . firmid: Matrix with firm fixed effect identifiers. Dimensions: N* x 1

    
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-                   
                        %---NON-MANDATORY INPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-
    
    . other_fe: The user can add more fixed effects than worker and firm
    ones but the connected set and the leave-one-out connected set are
    computed using id and firmid. The additional fixed effect identifiers
    are normalized to avoid multicolinearities.

    . mat_controls: Matrix of controls. Dimensions: N* x K. 
	
	. group: If the estimation of the second order moments and the variance 
	decomposition needs to be splitted into different groups the user needs 
	to provide this variable. The code accepts performing several groupings 
	G if a matrix is provided. Dimensions:  N* x G.	
    
    . n_lev: Number of simulations for leverage estimation. It is required when
    choosing type_hc equal to 'hc_2', 'hc_u', 'hc_u_match' or 'hc_u_clus'.
	A natural number. The default value is 300.

    . n_boot: Number of bootstrap simulations. This parameter governs the
    precision of the estimation. A natural number. The default value is 300.

    . type_hc: Can take values of 'hom', 'hc_0', 'hc_1', 'hc_2', 'hc_u', 
    'hc_u_match' or 'hc_u_clus' for different estimates of the covariance
    matrix. 'hom', 'hc_0', 'hc_1', 'hc_2', 'hc_u' are covariance matrix 
    estimators without match or cluster dependence. 
    'hc_u_match' assumes error dependence within the match and 'hc_u_clus'
    assumes error dependence within a user-defined cluster. The 
    The default value is 'hc_u_match'.

    . type_leave: Can take values of 'obs', 'worker' or
    'match'
	
	. cluster: vector defining the clustering of the errors. The 

    . ind_light: indicator to keep a light environment by clearing
    variables from the base workspace. The default is ind_light = 0.

    . mytol: Tolerance for pcg when solving the normal equations. The
    default value is 1e-6.

    . LdM_mom: it can take values of 1 or 0. If 1 then the algorithm also
    computes an extra vector taking the averages of worker fixed effects
    per firm. It then uses this extra vector to compute
    its variance and covariance with the other vectors. With this we could,
    for example, build a corrected correlation of worker fixed effects with
    the average of worker fixed effects within a firm as in Lopes de Melo
    (2018, JPE).
	
	. year: If LdM_mom is equal to 1, year is a required input.

    . filename: path and name to store the second order moment estimates and
    the variance decompositions. String. 
    If the user provides defines 'ind_export' as 1 with filename = 'example' and a grouping vector 'group', the output files will be: 
        example_plugin_estimates.csv 
        example_corrected_estimates.csv
        example_var_decomp_plugin.csv
        example_var_decomp_corrected.csv
        example_N_obs_group.csv

    . ind_export: indicator of printing csv files with results and number of
    observations per group if 'group' variable is provided. Can take values
    of 0 and 1. If 'filename' is not provided, the files will be printed in
    the current directory as: '_plugin_estimates.csv',
    '_corrected_estimates.csv', '_var_decomp_plugin.csv',
    _var_decomp_corrected.csv', and '_N_obs_group.csv' if 'group' is
    provided

    . v_filename_group: vector of grouping filenames to be added to main
    filename if group is a matrix

    
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%-             
    %% 					DESCRIPTION OF THE OUTPUTS
    %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- %-%-%- 
 
    . plugin: Plugin estimates of second order moments.  If the user provided a 
	'group' variable, the plugin estimates are also computed per group.

    . delta: Bootstrap estimate of the bias of second order moments. If the 
	user provided a 'group' variable, the bias is also estimated per group.

    . corrected: Bootstrap corrected second order moments. If the 
	user provided a 'group' variable, the correction is also estimated per group.

    . decomp_pi: Variance decomposition using the plugin estimates of
    second order moments. If the user provided a 'group' variable, the
    decomposition is also performed per group.

    . decomp_b: Variance decomposition using the bootstrap corrected
    estimates of second order moments. If the user provided a 'group' vector, the
    decomposition is also performed per group.

    . dimensions: Dimensions taking into account the normalizations

    . NT: # of person-year observations in the final sample where the
    second order moments are computed. NT can be below N* if the user did
    not provide a connected set or a leave-one-out connected set.

    . Stored CSVs: The function stores 4 csv files with the plugin
    estimates, the corrected estimates and the respective variance
    decompositions. The user can choose the starting filename for the
    output with the optional input 'filename'.

%}

    %% Initial variables

    % ----------------------------------------------------------------------
    %   START
    % ----------------------------------------------------------------------
    
    % Define default parameters
    default_param.controls = [];
    default_param.n_lev = 300;
    default_param.n_boot = 300;
    default_param.group = [];
    default_param.type_hc = 'hc_u_match';
    default_param.type_leave = []; % to be filled later on
    default_param.ind_light = 0;
    default_param.mytol = 1e-6;
    default_param.filename = '';
    default_param.v_filename_group = '';
    default_param.LdM_mom = 0;
    default_param.ind_export = 0;
    default_param.other_fe = [];

    % Assign default values if the variable is missing
    if nargin < 2
        error('More arguments needed');
    else
        if ~exist('mat_controls', 'var') || size(mat_controls,2)==0
            mat_controls = default_param.controls;
        end
        if ~exist('n_lev', 'var') || size(n_lev,2)==0
            n_lev = default_param.n_lev;
        end
        if ~exist('n_boot', 'var') || size(n_boot,2)==0
            n_boot = default_param.n_boot;
        end
        if ~exist('group', 'var') || size(group,2)==0
            group = default_param.group;
        end            
        if ~exist('type_hc', 'var') || size(type_hc,2)==0
            type_hc = default_param.type_hc;
        end           
        if ~exist('type_leave', 'var') || size(type_leave,2)==0
            type_leave = default_param.type_leave;
        end 
        if ~exist('ind_light', 'var') || size(ind_light,2)==0
            ind_light = default_param.ind_light;
        end   

         if ~exist('mytol', 'var') || size(mytol,2)==0
            mytol = default_param.mytol;
        end           
        if exist('LdM_mom', 'var')==0 || size(LdM_mom,2)==0
            LdM_mom = default_param.LdM_mom;
        end
        if exist('filename', 'var')==0 || size(filename,2)==0
            filename = default_param.filename;
        end
        if exist('v_filename_group', 'var')==0 || size(v_filename_group,2)==0
            v_filename_group = default_param.v_filename_group;
        end
        if exist('ind_export', 'var')==0 || size(ind_export,2)==0
            ind_export = default_param.ind_export;
        end
        if exist('other_fe', 'var')==0
            other_fe = default_param.other_fe;
        end
        if (exist('year', 'var')==0 || size(year,2)==0) && LdM_mom==1
            disp("A year vector is required to compute the average worker effect of coworkers in a year. The code cannot proceed.")
            error("Need a vector of year identifiers for LdM moment.")
        end

    end
    

    %% Generate cases
    if strcmp(type_hc,'hc_u_match')
        if (isempty(group) && LdM_mom == 0)
    
            disp("General options: Standard corrections with leave match out variance")
            
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT]= correction_match(y, id, firmid, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                    type_leave, ind_light, mytol, ind_export, filename);
    
        elseif (isempty(group) && LdM_mom == 1) 
    
            disp("General options: Standard correction with LdM moment")
                                                                                            
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_match_LdM(y, id, firmid, year, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                   type_leave, ind_light, mytol, ind_export, filename);

        elseif (~isempty(group) && LdM_mom == 0) 
          disp("General options: Corrections per multi-group with leave match out variance")
            
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_match_multi_group(y, id, firmid, group, other_fe, mat_controls, n_lev, n_boot,...
                type_leave, ind_light, mytol, ind_export, filename, v_filename_group);
             
        elseif (~isempty(group) && LdM_mom == 1) 
            disp("General options: Corrections per multi-group with LdM moment")
    
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_match_multi_group_LdM(y, id, firmid, year, group, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                    type_leave, ind_light, mytol, ind_export, filename, v_filename_group);

        end

    elseif strcmp(type_hc,'hc_u_clus')
        if (isempty(group) && LdM_mom == 0)
    
            disp("General options: Corrections with leave cluster out variance")
             [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_cluster(y, id, firmid, cluster, other_fe, mat_controls, n_lev,...    
                                                                                                n_boot, ind_light, mytol, ind_export, filename);

        elseif (isempty(group) && LdM_mom == 1)
        
             disp("General options: Corrections with leave cluster out variance with LdM moment")
             [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_cluster_LdM(y, id, firmid, cluster, year, other_fe, mat_controls, n_lev,...    
                                                                                                n_boot, ind_light, mytol, ind_export, filename);
          
        elseif (~isempty(group) && LdM_mom == 0) 
    
            disp("General options: Corrections per group with leave cluster out variance")
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_cluster_multi_group(y, id, firmid, cluster, group, other_fe, mat_controls, n_lev,...    
                                                                                                n_boot, ind_light, mytol, ind_export, filename, v_filename_group);
        
        elseif (~isempty(group) && LdM_mom == 1) 
    
            disp("General options: Corrections per group with leave cluster out variance with LdM moment")
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_cluster_multi_group_LdM(y, id, firmid, cluster, year, group, other_fe, mat_controls, n_lev,...    
                                                                                                n_boot, ind_light, mytol, ind_export, filename, v_filename_group);
      
        end

    else   % diagonal covariance matrix
        if (isempty(group) && LdM_mom == 0)
      
            disp("General options: Standard correction without group decompositions")

            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction(y, id, firmid, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                    type_hc, type_leave, ind_light, mytol, ind_export, filename);   

        elseif (isempty(group) && LdM_mom == 1)
           
            disp("General options: Standard correction with LdM moment")

            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_LdM(y, id, firmid, year, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                    type_hc, type_leave, ind_light, mytol, ind_export, filename);
    
        elseif (~isempty(group) && LdM_mom == 0) 
    
            disp("General options: Corrections per group")
            
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_multi_group(y, id, firmid, group, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                    type_hc, type_leave, ind_light, mytol, ind_export, filename, v_filename_group);
      
        elseif (~isempty(group) && LdM_mom == 1) 
    
            disp("General options: Corrections per group with LdM moment")
    
            [plugin, delta, corrected, decomp_pi, decomp_b, dimensions, NT] = correction_multi_group_LdM(y, id, firmid, year, group, other_fe, mat_controls, n_lev, n_boot,...
                                                                                                    type_hc, type_leave, ind_light, mytol, ind_export, filename, v_filename_group);

    
        end

    end

