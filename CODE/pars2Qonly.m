function [kinfQ,K1Q_X] = pars2Qonly(pars)
% This function distributes Q parameters from a stacked vector

% distribute parameters
kinfQ = pars(1);
K1Q_X = diag(pars(2:end));

end