function [hc_adj,n_problems] = diagnostic_par(correction_lev, NT, X, xx, L, mii, mytol,lap)
    
    %% Diagnosis of problems in the leverage estimation
           
    % Check hc after non-linearity correction
    hc_check = 1./mii.*correction_lev;
    
    % Problem 1: leverage below lower bound (hii < 1/NT)
    % Equivalent for 1/mii < NT/(NT-1)
    diagnos_min = (hc_check<(NT./(NT-1))).*(hc_check>=0);
    
    % Problem 2: leverage larger than upper bound (hii > 1)
    % Equivalent to 1/mii<0 
    diagnos_max = (hc_check<0);
    
    %%% Overall problem indicator
    problems = ((diagnos_min>0) | (diagnos_max>0)); % indicator of problems
    
    n_problems = sum(problems);
    
    disp(['Problems in leverage estimation: ',num2str(n_problems)])
    
    %% Direct computation of problematics
    
    if n_problems>0
        disp('---- Origin of problems -----')
        disp(['  Below minimum: ',num2str(sum(diagnos_min))])
        disp(['  Above maximum: ',num2str(sum(diagnos_max))])
        disp('Direct computation of problematics')
    
    pos_problem = find(problems);
    
    lev_p = zeros(n_problems,1);
    
    % Pass variables to workers
    X_par = parallel.pool.Constant(X);
    xx_par = parallel.pool.Constant(xx);
    pos_problem_par = parallel.pool.Constant(pos_problem);


        if lap==1
        
            L_par = parallel.pool.Constant(@() L);
                        
            parfor j=1:n_problems
        
                v_prime = zeros(NT,1);
                v_prime(pos_problem_par.Value(j)) = 1;
                        
                % Run the regression
                [b_lev,~]=pcg(xx_par.Value,X_par.Value'*v_prime,mytol,300,L_par.Value);
                
                % Fitted value
                xb0 = X_par.Value*b_lev;
                lev_p(j) = xb0(pos_problem_par.Value(j));     % fitted values of the problematic ones
            
            end
            
            % assign real leverage correction to problematics
            hc_adj = hc_check;
            hc_adj(problems) = 1./(1-lev_p);
            
        else
            
            L_par = parallel.pool.Constant(L);
                        
            parfor j=1:n_problems
   
                v_prime = zeros(NT,1);
                v_prime(pos_problem_par.Value(j)) = 1;
                        
                % Run the regression
                [b_lev,~]=pcg(xx_par.Value,X_par.Value'*v_prime,mytol,300,L_par.Value,L_par.Value');
                
                % Fitted value
                xb0 = X_par.Value*b_lev;
                lev_p(j) = xb0(pos_problem_par.Value(j));     % fitted values of the problematic ones
       
            end
            
            % assign real leverage correction to problematics
            hc_adj = hc_check;
            hc_adj(problems) = 1./(1-lev_p);           
            
        end
        
        disp('Done!')
        
    else
    
        hc_adj = hc_check;
      
    end
        
end