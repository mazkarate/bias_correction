function [Pii_bar, Mii_bar, correction_lev] = est_lev(n_lev, NT_mover, X, xx, L , id, mytol)

Pii_mover = 0;
Mii_mover = 0;
Mii_sq_mover = 0;
Pii_Mii_mover = 0;

parfor j=1:n_lev

    y_lev = randsample([-1;1], NT_mover, true);
        
    % Estimate regression
    [b_lev,~]=pcg(xx, X'*y_lev, mytol, 300, L);
                
    y_hat=X*b_lev;     % fitted values
        
    %Collect estimates for leverage (Pii) and 1 - leverage (Mii)
    Pii_mover = Pii_mover + y_hat.^2;
            
    Mii_mover = Mii_mover + ((y_lev - y_hat).^2);
    Mii_sq_mover = Mii_sq_mover + ((y_lev - y_hat).^4);
        
    Pii_Mii_mover = Pii_Mii_mover + ((y_hat.^2).*((y_lev-y_hat).^2));
        
end
%Get averages
Pii_mover = Pii_mover./n_lev;
   
Mii_mover = Mii_mover./n_lev;
Mii_sq_mover = Mii_sq_mover./n_lev;
        
Pii_Mii_mover = Pii_Mii_mover./n_lev;

%% Rescale and bias correction of leverage estimates

%Here we rescale the estimates of both Pii and Mii which guarantees that
%both lie within [0,1]. 

%We do our own bias correction for the non-linear estimate of 1/M_ii. Check
%the online appendix for details.

Pii_bar_mover             = Pii_mover./(Pii_mover+Mii_mover); 	
Mii_bar_mover             = 1-Pii_bar_mover;
var_m               = (1./(n_lev-1)).*(Mii_sq_mover-Mii_mover.^2);
cov_hm              = (1./(n_lev-1)).*(Pii_Mii_mover - Pii_mover.*Mii_mover);
correction_lev 		= 1 - var_m.*Pii_bar_mover./Mii_bar_mover.^2 + cov_hm./Mii_bar_mover;

%% Estimates of movers and stayers
T = accumarray(id,1);
T = T(id); % NT size vector

T_stayer = T(NT_mover + 1: end);
P_ii_stayer = 1./T_stayer;
M_ii_stayer = 1-P_ii_stayer;

%%% estimates for all
Pii_bar = [Pii_bar_mover; P_ii_stayer];
Mii_bar = [Mii_bar_mover; M_ii_stayer];


end

