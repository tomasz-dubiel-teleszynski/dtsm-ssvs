function K1Q_X = mlgitinvs(K1Q_X)
% This function transforms back diagonal elements of K1Q_X using inverse
% logit which preseves their decreasing order

% get diagonal
diagonal = diag(K1Q_X);

% rescaling factor 
rescaler = 1;

% rescale
diagonal = diagonal/rescaler;

% transform back
diagonal = - log( - diagonal);
 
K1Q_X = diag(cumsum(diagonal));

end
