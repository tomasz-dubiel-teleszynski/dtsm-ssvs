function [Llk,LlkP,LlkQ] = getLlk1(Y,W,cP,mats,dt,...
                                   kinfQ,...
                                   K1Q_X,...
                                   lam0,...
                                   lam1,...
                                   Sigma_cP,...
                                   sige)
% This function calculates model log likelihood as sum of Q and P log
% likelihoods for 1 observation

% get Q parameters for observed factors cP and Ricatti loadings for
% observed yields from factor loading matrix and Q parameters for
% unobserved factors X
[K0Q_cP,K1Q_cP,AcP,BcP] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);

% get Q log likelihood for 1 observation
LlkQ = getLlkQ1(Y,W,cP,AcP,BcP,sige);

% get P parameters for observed factors cP from their Q parameters
[K0P_cP,I_K1P_cP] = getPdyn(K0Q_cP,K1Q_cP,lam0,lam1);

% get P log likelihood
LlkP = getLlkP(cP,K0P_cP,I_K1P_cP,Sigma_cP);

% get model log likelihood as sum of P and Q log likelihoods
Llk = LlkP + LlkQ;

end

