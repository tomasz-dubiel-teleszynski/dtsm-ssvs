function logPrior = uvnlogPrior(pars,EX,VAR)
% This function calculates log prior density as sum of univariate log N(EX,VAR)
% densities up to proportionality for vector of parameters in pars

% sum of univariate log N(EX,VAR) densities up to proportionality
logPrior = sum(-0.5*log(VAR) - ((pars - EX).^2)./(2*VAR));

end