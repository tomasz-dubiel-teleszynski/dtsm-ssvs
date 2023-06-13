function pars = Qonly2pars(kinfQ,K1Q_X)
% This function stacks Q parameters into a vector

% stack Q parameters
pars = [kinfQ;diag(K1Q_X)];

end