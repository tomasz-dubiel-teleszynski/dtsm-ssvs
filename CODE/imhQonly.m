function [kinfQ,K1Q_X,logPost] = imhQonly(logPost,...
                                          kinfQ,...
                                          K1Q_X,...
                                          parsQonly_mle,...
                                          covQonly_mle,...
                                          EX_Qonly,...
                                          VAR_Qonly,...
                                          EX_Sigma_cP,...
                                          VAR_Sigma_cP,...
                                          Y,W,cP,mats,dt,lam0,lam1,Sigma_cP,sige,...
                                          varargin)
% This function implements independence Metropolis-Hastings sampler jointly
% for kinfQ and K1Q_X assuming univariate N(EX_Qonly,VAR_Qonly) prior for
% each Q parameter
                              
% scale and transform 
parsQonly = Qonly2pars(scalekinfQ(kinfQ),mlgits(K1Q_X));

% propose new kinfQ and K1Q_X
nu = 5;
pars1Qonly = MVT_RND(parsQonly_mle,(nu - 2)/nu*covQonly_mle,nu,1);

% rescale and transform back
[kinfQ1,K1Q_X1] = pars2Qonly(pars1Qonly);
kinfQ1 = rescalekinfQ(kinfQ1);
K1Q_X1 = mlgitinvs(K1Q_X1);

% singularity check of matrix WBX in "jszloadings.m" after proposal
if isWBXsingular(W,mats,dt,kinfQ1,K1Q_X1)
   return;
end

% tempering
if isempty(varargin)
    fi = 1;
else
    fi = varargin{1};
end

% log posterior density
% at kinfQ1 and K1Q_X1
logPost1 = getlogPostQonly(pars1Qonly,...
                           EX_Qonly,...
                           VAR_Qonly,...
                           EX_Sigma_cP,...
                           VAR_Sigma_cP,...
                           Y,W,cP,mats,dt,kinfQ1,K1Q_X1,lam0,lam1,Sigma_cP,sige,...
                           fi);

% at kinfQ and K1Q_X - THIS CAN BE COMMENTED OUT TO SAVE TIME
% logPost = getlogPostQonly(parsQonly,...
%                           EX_Qonly,...
%                           VAR_Qonly,...
%                           EX_Sigma_cP,...
%                           VAR_Sigma_cP,...
%                           Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige,...
%                           fi);
                                        
% log proposal density
% at pars1Qonly
logProp1 = mvtlogProp(pars1Qonly,parsQonly_mle,covQonly_mle,nu);

% at parsQonly
logProp = mvtlogProp(parsQonly,parsQonly_mle,covQonly_mle,nu);

% log acceptance probability
logAcceptProb = logPost1 - logPost + logProp - logProp1;

% check if transition matrix under P for observed factors is not
% explosive that is check if all its eigenvalues are below 1 in absolute
% value
[~,K1Q_cP] = jszloadings(W,mats,dt,kinfQ1,K1Q_X1,Sigma_cP);
IN = eye(numel(lam0));
tmp = eig(K1Q_cP + IN + lam1);
% check if conditions for eigenvalues hold [if they lie in (-1;0)] and
% if they are sorted
tmp1 = diag(K1Q_X1);
% independence Metropolis accept-reject step
if (log(rand()) < logAcceptProb) &&...
        all(abs(tmp) < 1) &&...
        all(tmp1 < -eps) &&...
        all(tmp1 > -1 + eps) &&...
        all(diff(tmp1) < -eps)
    kinfQ = kinfQ1;
    K1Q_X = K1Q_X1;
    logPost = logPost1;
end

end