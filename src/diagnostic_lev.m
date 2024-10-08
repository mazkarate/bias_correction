function [hc_adj, n_problems] = diagnostic_lev(correction_lev, NT, X, xx, L, Mii, mytol)
    
    %% Diagnosis of problems in the leverage estimation
           
    % Check hc after non-linearity correction
    hc_check = 1./Mii.*correction_lev;
    
    % Problem 1: leverage below lower bound (Pii < 1/NT)
    % Equivalent for 1/Mii < NT/(NT-1)
    diagnos_min = (hc_check<(NT./(NT-1))).*(hc_check>=0);
    
    % Problem 2: leverage larger than upper bound (Pii > 1)
    % Equivalent to 1/Mii<0 
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
        pos_problem_par = parallel.pool.Constant(pos_problem);

        parfor j=1:n_problems
        
            v_prime = zeros(NT,1);
            v_prime(pos_problem_par.Value(j)) = 1;
                        
            % Run the regression
            [b_lev,~]=pcg(xx,X'*v_prime,mytol,300,L);
                
            % Fitted value
            xb0 = X*b_lev;
            lev_p(j) = xb0(pos_problem_par.Value(j));     % fitted values of the problematic ones
            
        end
            
        % assign real leverage correction to problematics
        hc_adj = hc_check;
        hc_adj(problems) = 1./(1-lev_p);
               
        disp('Done!')
        
    else
    
        hc_adj = hc_check;
      
    end
        
end