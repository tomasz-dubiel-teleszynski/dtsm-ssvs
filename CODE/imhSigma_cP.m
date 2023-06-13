function [Sigma_cP,logPost] = imhSigma_cP(Sigma_cP,...
                                          parsSigma_cP_mle,...
                                          covSigma_cP_mle,...
                                          EX_Sigma_cP,...
                                          VAR_Sigma_cP,...
                                          EX_Qonly,...
                                          VAR_Qonly,...
                                          Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,sige,...
                                          varargin)
% This function implements independence Metropolis-Hastings sampler for
% Sigma_cP assuming univariate N(EX_Sigma_cP,VAR_Sigma_cP) prior for each
% of its non-zero elements

% transform and scale                 
parsSigma_cP = Sigma_cP2pars(logd(Sigma_cP));

% propose new Sigma_cP
nu = 5;
pars1Sigma_cP = MVT_RND(parsSigma_cP_mle,(nu - 2)/nu*covSigma_cP_mle,nu,1);

% transform back and rescale
Sigma_cP1 = expd(pars2Sigma_cP(pars1Sigma_cP));

% tempering
if isempty(varargin)
    fi = 1;
else
    fi = varargin{1};
end

% log posterior density
% at Sigma_cP1
logPost1 = getlogPostSigma_cP(pars1Sigma_cP,...
                              EX_Sigma_cP,...
                              VAR_Sigma_cP,...
                              EX_Qonly,...
                              VAR_Qonly,...
                              Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,Sigma_cP1,sige,...
                              fi);
                      
% at Sigma_cP
logPost = getlogPostSigma_cP(parsSigma_cP,...
                             EX_Sigma_cP,...
                             VAR_Sigma_cP,...
                             EX_Qonly,...
                             VAR_Qonly,...
                             Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige,...
                             fi);
                 
% log proposal density
% at pars1Sigma_cP
logProp1 = mvtlogProp(pars1Sigma_cP,parsSigma_cP_mle,covSigma_cP_mle,nu);

% at parsSigma_cP
logProp  = mvtlogProp(parsSigma_cP,parsSigma_cP_mle,covSigma_cP_mle,nu);

% log acceptance probability
logAcceptProb = logPost1 - logPost + logProp - logProp1;

% check if transition matrix under P for observed factors is not
% explosive that is check if all its eigenvalues are below 1 in absolute
% value
[~,K1Q_cP] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP1);
IN = eye(numel(lam0));
tmp = eig(K1Q_cP + IN + lam1);
% independence Metropolis-Hastings proposal accept-reject step
if (log(rand()) < logAcceptProb) &&...
        all(abs(tmp) < 1)
    Sigma_cP = Sigma_cP1;
    logPost = logPost1;
end

end