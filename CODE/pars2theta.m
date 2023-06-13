function [K1Q_X,Sigma_cP] = pars2theta(pars)

a = numel(pars);
N = (sqrt(9 + 8*a) - 3)/2;

K1Q_X = diag(pars(1:N));
Sigma_cP = pars2Sigma_cP(pars(N+1:end));

end