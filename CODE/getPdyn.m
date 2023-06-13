function [K0P_cP,I_K1P_cP] = getPdyn(K0Q_cP,K1Q_cP,lam0,lam1)
% This function calculates P parameters of VAR(1) for observed factors cP
% from Q parameters of VAR(1) for these factors

K0P_cP = K0Q_cP + lam0;
K1P_cP = K1Q_cP + lam1;
I_K1P_cP = eye(numel(lam0)) + K1P_cP;

end