function hessLlk = hessLlkQonly(pars,Y,W,cP,mats,dt,lam0,lam1,Sigma_cP,sige)
% This function calculates model log likelihood at Q parameters

% distribute stacked Q parameters
[kinfQ,K1Q_X] = pars2Qonly(pars);
% rescale
kinfQ = rescalekinfQ(kinfQ);
% transform back
K1Q_X = mlgitinvs(K1Q_X);

% calculate model log likelihood
hessLlk = getLlk(Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige);
hessLlk = - hessLlk;

end

