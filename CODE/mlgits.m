function K1Q_X = mlgits(K1Q_X)
% This function performs a logit transformation on diagonal elements of
% K1Q_X and preserves their decreasing order

% get diagonal
diagonal = diag(K1Q_X);

% scaling factor
scaler = 1;

% scale
diagonal = diagonal*scaler; 

% transform
K1Q_X = diag( - exp( - [diagonal(1);diff(diagonal)]));

end