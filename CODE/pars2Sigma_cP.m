function Sigma_cP = pars2Sigma_cP(pars)
% This function distributes Sigma_cP non-zero elements from a stacked
% vector

% distribute parameters
a = numel(pars);
N = (sqrt(1 + 8*a) - 1)/2; 
if N == 3
    Sigma_cP =  [pars(1),      0,      0;
                 pars(2),pars(4),      0;
                 pars(3),pars(5),pars(6)];
elseif N == 2
    Sigma_cP =  [pars(1),      0;
                 pars(2),pars(3)]; 
end
         
end