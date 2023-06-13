function Sigma_cP = expdmle(Sigma_cP)
% This function exponentiates diagonal elements of Sigma_cP and rescales
% its other non-zero elements

% exponentiate diagonal elements
Sigma_cP(eye(size(Sigma_cP))~=0) = exp(diag(Sigma_cP));

% rescaling factor
rescaler = 1e6;

% rescale
N = numel(diag(Sigma_cP));
if N == 3
    Sigma_cP(2:3,1) = Sigma_cP(2:3,1)/rescaler;
    Sigma_cP(3,2) = Sigma_cP(3,2)/rescaler;
elseif N == 2
    Sigma_cP(2,1) = Sigma_cP(2,1)/rescaler;
end

end

