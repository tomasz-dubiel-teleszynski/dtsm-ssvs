function [W,cP] = getWcP(Y,N)
% This function calculates PCA with normalization as in Bauer (2016)

% do PCA
[U,~,~] = svd(cov(Y));  
J = size(Y,2);
W = U(:,1:N)';

% do normalizations
tmp = sum(W(1,:));
W(1,:) =  W(1,:)/tmp; % normalize level so that loadings sum to one
tmp = W(2,J) - W(2,1);
W(2,:) = W(2,:)/tmp; % normalize slope so that difference between loadings on long and short yield is one
if N > 2
    % normalize curvature
    tmp = W(3,J) + W(3,1) - 2*W(3,floor(J/2));
    W(3,:) = W(3,:)/tmp;    
end

% extract observed factors from yields
cP = Y*W';

end

