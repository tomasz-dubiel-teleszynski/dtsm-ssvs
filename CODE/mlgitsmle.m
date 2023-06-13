function K1Q_X = mlgitsmle(K1Q_X)
% This function performs a logit transformation on diagonal elements of
% K1Q_X and preserves their decreasing order

% get diagonal
diagonal = diag(K1Q_X);

% transform
N = numel(diagonal);
if N == 3
    K1Q_X(eye(size(K1Q_X))~=0) = [logit(diagonal(1),-1,0),...
                                  logit(diagonal(2),-1,diagonal(1)),...
                                  logit(diagonal(3),-1,diagonal(2))];
elseif N == 2
    K1Q_X(eye(size(K1Q_X))~=0) = [logit(diagonal(1),-1,0),...
                                  logit(diagonal(2),-1,diagonal(1))];
end

end