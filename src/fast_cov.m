% Betthauser, J., 2017, Email: jbettha1@jhu.edu
% PhD candidate, Electrical and Computer Engineering
% Johns Hopkins University, Brain-Computer Interfaces and Neuroprostheses Lab
%% fast_cov(): Covariance calculation 15-35% faster than MATLAB built-in cov()
%       Input X : M x N ( # features x # samples)
%       Output Xcov : M x M covariance matrix
% NOTE: Confirmed that output of fast_cov(X) == cov(X')
function [ Xcov ] = fast_cov( X )
    mu = mean(X,2); len = size(X,2);
    Xcentered = X - repmat(mu,1,len);
    Xcov = Xcentered * Xcentered' ./ (len - 1);
end