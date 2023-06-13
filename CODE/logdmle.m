function Sigma_cP = logdmle(Sigma_cP)
% This function log transforms diagonal elements of Sigma_cP and scales its
% other non-zero elements

% diagonal log transform
Sigma_cP(eye(size(Sigma_cP))~=0) = log(diag(Sigma_cP)); 

% scaling factor
scaler = 1e6;

% scale
N = numel(diag(Sigma_cP));
if N == 3   
    Sigma_cP(2:3,1) = scaler*Sigma_cP(2:3,1);
    Sigma_cP(3,2) = scaler*Sigma_cP(3,2);
elseif N == 2
    Sigma_cP(2,1) = scaler*Sigma_cP(2,1);
end

end

