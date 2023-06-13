function pars = Sigma_cP2pars(Sigma_cP)
% This function stacks non-zero elements of Sigma_cP into a vector

% stacking columnwise, omitting above diagonal elements
N = numel(diag(Sigma_cP));
if N == 3
    pars = [Sigma_cP(1:3,1);Sigma_cP(2:3,2);Sigma_cP(3,3)];
elseif N == 2
     pars = [Sigma_cP(1:2,1);Sigma_cP(2,2)];
end
    
end