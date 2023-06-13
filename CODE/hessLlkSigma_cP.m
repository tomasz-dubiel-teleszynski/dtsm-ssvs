function hessLlk = hessLlkSigma_cP(pars,Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,sige)
% This function calculates model log likelihood at Sigma_cP

% distribute stacked elements of Sigma_cP, epxonentiate diagonal and
% rescale its other non-zero elements
Sigma_cP = expd(pars2Sigma_cP(pars));

% calculate model log likelihood
hessLlk = getLlk(Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige);
hessLlk = - hessLlk;
    
end

