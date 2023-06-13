function covSigma_cP = getcovSigma_cPi(eps,lift,Sigma_cP,Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,sige)
% This function calculates covariance matrix for Sigma_cP at the ML
% estimate as inverse of negative hessian (Fisher information matrix) and
% increments its diaognal elements if specified

% stack non-zero elements of Sigma_cP into a vector, transform and scale
parsSigma_cP = Sigma_cP2pars(logd(Sigma_cP));

% get hessian at Sigma_cP
hessSigma_cP = numhessi(@hessLlkSigma_cP,parsSigma_cP,eps,Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,sige);

% calculate covariance matrix as minus inverse of the Hessian
covSigma_cP = makepd(hessSigma_cP\eye(numel(parsSigma_cP)));

% increment diagonal elements
idx = (eye(numel(parsSigma_cP)) == 1);
covSigma_cP(idx) = (1 + lift)*covSigma_cP(idx);
covSigma_cP = makepd(covSigma_cP);

end