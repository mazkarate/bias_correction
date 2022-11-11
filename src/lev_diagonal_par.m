function [hii_bar, mii_bar, correction_lev] = lev_diagonal_par(n_lev, NT, X, xx, L , mytol, lap)
        
hii = 0;
%hii_sq = 0;
mii = 0;
mii_sq = 0;
hii_mii = 0;

if lap==1

    % Pass variables to workers
    X_par = parallel.pool.Constant(X);
    xx_par = parallel.pool.Constant(xx);
    L_par = parallel.pool.Constant(@() L);
    
    % Pass variables to workers 
    parfor j=1:n_lev
        
        y_lev = randsample([-1;1],NT,true);
        
        % Estimate regression
        [b_lev,~]=pcg(xx_par.Value,X_par.Value'*y_lev,mytol,300,L_par.Value);
        
        % Fitted value to the square
        y_hat=X_par.Value*b_lev;     % fitted values
        
        %Collect estimates for leverage (hii) and 1 - leverage (mii)
        hii = hii + y_hat.^2;
        %hii_sq = hii_sq + y_hat.^4;
        
        mii = mii + ((y_lev - y_hat).^2);
        mii_sq = mii_sq + ((y_lev - y_hat).^4);
        
        hii_mii = hii_mii + ((y_hat.^2).*((y_lev-y_hat).^2));
        
    end

else

    % Pass variables to workers
    X_par = parallel.pool.Constant(X);
    xx_par = parallel.pool.Constant(xx);
    L_par = parallel.pool.Constant(L);
    
    % Pass variables to workers 
    parfor j=1:n_lev
        
        y_lev = randsample([-1;1],NT,true);
        
        % Estimate regression
        [b_lev,~]=pcg(xx_par.Value,X_par.Value'*y_lev,mytol,300,L_par.Value,L_par.Value');

        % Fitted value to the square
        y_hat=X_par.Value*b_lev;     % fitted values
        
        %Collect estimates for leverage (hii) and 1 - leverage (mii)
        hii = hii + y_hat.^2;
        %hii_sq = hii_sq + y_hat.^4;
        
        mii = mii + ((y_lev - y_hat).^2);
        mii_sq = mii_sq + ((y_lev - y_hat).^4);
        
        hii_mii = hii_mii + ((y_hat.^2).*((y_lev-y_hat).^2));
        
    end

end

%Get averages
hii = hii./n_lev;
%hii_sq = hii_sq./n_lev;
   
mii = mii./n_lev;
mii_sq = mii_sq./n_lev;
        
hii_mii = hii_mii./n_lev;

%% Rescale and bias correction of leverage estimates

%Here we rescale the estimates of both hii and mii which guarantees that
%both lie within [0,1]. We then do a bias correction for the non-linear
%estimation of hii and mii.

%We follow the improved code for this leverage estimation by KSS.

%Reference:
%Kline, Patrick, Raffaele Saggio, and Mikkel Sølvsten (2021). Improved 
%stochastic approximation of regression leverages for bias correction of 
%variance components. 
%https://www.dropbox.com/s/i28yvzae2tnp2tl/improved_JLA.pdf?dl=1

%We do our own bias correction for the non-linear estimate of 1/m_ii. Check
%the online appendix for details.

hii_bar             = hii./(hii+mii); 	
mii_bar             = 1-hii_bar;
var_m               = (1./(n_lev-1)).*(mii_sq-mii.^2);
cov_hm              = (1./(n_lev-1)).*(hii_mii - hii.*mii);
correction_lev 		= 1 - var_m.*hii_bar./mii_bar.^2 + cov_hm./mii_bar;

end

