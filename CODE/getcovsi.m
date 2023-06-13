function [covSigma_cP_mle,covQonly_mle] = getcovsi(kinfQ_mle,...
                                                   K1Q_X_mle,...
                                                   lam0_mle,...
                                                   lam1_mle,...
                                                   Sigma_cP_mle,...
                                                   sige_mle,...
                                                   Y,W,cP,mats,dt)
% This function calculates covariance matrices for Sigma_cP and for Q
% paramters at the ML estimates as inverse of negative hessians
% (Fisher information matrices) and increments their diagonal elements
% as specified below
     
% transform and scale Sigma_cP_mle for independence Metropolis-Hastings
% sampler
parsSigma_cP_mle = Sigma_cP2pars(logd(Sigma_cP_mle));

% transform and scale kinfQ_mle and K1Q_X_mle for independence
% Metropolis-Hastings sampler
parsQonly_mle = Qonly2pars(scalekinfQ(kinfQ_mle),mlgits(K1Q_X_mle));

%% parameter hessians and covariance matrices at their MLEs
% covariance matrix for ML estimate of Sigma_cP
% tolerance
epsSigma_cP = 10.^(real(floor(log10(parsSigma_cP_mle))) - 4);

% diagonal increment
liftSigma_cP = 0.0;

% covariance matrix
covSigma_cP_mle = getcovSigma_cPi(epsSigma_cP,...
                                  liftSigma_cP,...
                                  Sigma_cP_mle,...
                                  Y,W,cP,mats,dt,...
                                  kinfQ_mle,...
                                  K1Q_X_mle,...
                                  lam0_mle,...
                                  lam1_mle,...
                                  sige_mle);

% covariance matrix for ML estimates of kinfQ and K1Q_X
% tolerance
epsQonly = 10.^(real(floor(log10(parsQonly_mle))) - 4);

% diagonal increment
liftQonly = 0.0;
                              
% covariance matrix
covQonly_mle = getcovQonlyi(epsQonly,...
                            liftQonly,...
                            kinfQ_mle,...
                            K1Q_X_mle,...
                            Y,W,cP,mats,dt,...
                            lam0_mle,...
                            lam1_mle,...
                            Sigma_cP_mle,...
                            sige_mle);
                        
end