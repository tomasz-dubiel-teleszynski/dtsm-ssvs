function LlkP = getLlkP(cP,K0P_cP,I_K1P_cP,Sigma_cP)
% This function calculates P log likelihood

% get number of observations and factors
[T,N]= size(cP);

% calculate P log likelihood
e = cP(2:T,:) - (repmat(K0P_cP',T-1,1) + cP(1:T-1,:)*I_K1P_cP');
Omega_cP = Sigma_cP*Sigma_cP';
llkP = - 0.5*N*log(2*pi) - 0.5*log(det(Omega_cP)) - 0.5*diag((e/Omega_cP)*e'); 
LlkP = sum(llkP);

end