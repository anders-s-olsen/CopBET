function [Ct, Ht] = DCC(dat,iter)
% function [Ct, Ht] = DCC(dat)
%
% Estimate a multivariate GARCH model using the DCC estimator of Engle and Sheppard
% 
% INPUTS:
%
%      dat          Zero mean T by p matrix
% 
% OUTPUTS:
%
%      Ct           p by p by T array of conditional correlations
%      Ht           p by p by T array of conditional covariances
%
%
% File created by Martin Lindquist on 07/22/14
% Last update: 07/22/14

if isempty(iter)
    iter=1000;
end

[T,p] = size(dat);

[~, ~, Ht] = dcc_mvgarch(dat,1,1,1,1,iter);        % Ht is the dyanamic covariance matrix

% Compute dynamic correlation matrix Ct
Ct = Ht;
for i=1:T,
    Ct(:,:,i) = Ht(:,:,i)./sqrt(diag(Ht(:,:,i))*diag(Ht(:,:,i))');
end
