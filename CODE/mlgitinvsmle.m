function K1Q_X = mlgitinvsmle(K1Q_X)
% This function transforms back diagonal elements of K1Q_X using inverse
% logit which preseves their decreasing order

% get diagonal
diagonal = diag(K1Q_X);

% transform back
N = numel(diagonal);
if N == 3
    phi1 = invlogit(diagonal(1),-1,0);
    phi2 = invlogit(diagonal(2),-1,phi1);
    phi3 = invlogit(diagonal(3),-1,phi2);
    K1Q_X(eye(size(K1Q_X))~=0) = [phi1,phi2,phi3];
elseif N == 2
    phi1 = invlogit(diagonal(1),-1,0);
    phi2 = invlogit(diagonal(2),-1,phi1);
    K1Q_X(eye(size(K1Q_X))~=0) = [phi1,phi2];
end

end
