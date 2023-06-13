function olsSigma_cP = getolsSigma_cP(cP)
% This function calculates OLS/ML estimate of Sigma_cP from observed
% factors cP

T = size(cP,1);
cP1 = cP(1:end-1,:);
const = ones(T-1,1);
Z = [const, cP1];
Y = cP(2:end,:);
bhat = (Z'*Z)\(Z'*Y);
errors = Y - Z*bhat;
olsOmega_cP = 1/(T-2)*(errors')*errors; 
olsSigma_cP = chol(olsOmega_cP,'lower');

end