function [out,acr_,out_] = mhwg2ibisiSSVS(Niter,...
                                          in,...
                                          dims,...
                                          parsSigma_cP_mle,...
                                          parsQonly_mle,...
                                          covSigma_cP_mle,...
                                          covQonly_mle,...
                                          alpha_sige2,...
                                          beta_sige2,...
                                          EX_lam0lam1,...
                                          EX_Sigma_cP,...
                                          EX_Qonly,...
                                          VAR_lam0lam1,...
                                          VAR_Sigma_cP,...
                                          VAR_Qonly,...
                                          gamma,...
                                          Y,W,cP,mats,dt,...
                                          varargin)
% This function implements MCMC jittering step for IBIS algorithm

%%  initialization at current particle
[kinfQ,K1Q_X,lam0,lam1,Sigma_cP] = pars2svms(in,dims);

%% miscellaneous
% get number of factors
N = numel(diag(K1Q_X));
% identity matrix
IN = eye(N);
% get number of maturities
J = numel(mats);
% to annualize maturities
n_per = 12;
% # of large Ricatti loadings
JL = mats(end)/n_per;

% tempering
if isempty(varargin)
    fi = 1;
else
    if numel(varargin) == 2
        fi = varargin{2};
    else
        fi = 1;
    end
end

%% containers for posterior draws
out_eig_I_K1Q_cP = zeros(N,Niter);
out_eig_I_K1P_cP = zeros(N,Niter);
out_AcP          = zeros(J,Niter);
out_BcP          = zeros(N*J,Niter);
out_AX           = zeros(J,Niter);
out_BX           = zeros(N*J,Niter);
out_AXL          = zeros(JL,Niter);
out_BXL          = zeros(N*JL,Niter);

out_kinfQ    = zeros(1,Niter);
out_K1Q_X    = zeros(N,Niter);
out_lam0     = zeros(N,Niter);
out_lam1     = zeros(N*N,Niter);
out_Sigma_cP = zeros((N*N + N)/2,Niter);
out_sige     = zeros(1,Niter);

out_gamma      = zeros(N*N + N,Niter);
out_gamma_prob = zeros(N*N + N,Niter);

%% container for log posterior
out_logPost = zeros(1,Niter);

%% MCMC sampling
for iter = 1:Niter    
    % ### Gibbs step for sige
    sige = sqrt(sige2gs(alpha_sige2,...
                        beta_sige2,...
                        Y,W,cP,mats,dt,kinfQ,K1Q_X,Sigma_cP,...
                        fi));
    
    % ### Gibbs step for lam0 and lam1 with non-zero restrictions
    [lam0,lam1] = lam0lam1gsrSSVS(lam0,lam1,...
                                  gamma,...
                                  EX_lam0lam1,...
                                  VAR_lam0lam1,...
                                  W,cP,mats,dt,kinfQ,K1Q_X,Sigma_cP,...
                                  fi);
                              
    % ### Gibbs step for gamma with Bernoulli posterior
    [gamma, prob] = gammagsSSVS(lam0,lam1,...
                                gamma,...
                                EX_lam0lam1,...
                                VAR_lam0lam1);
                              
    % ### independence Metropolis-Hastings step for Sigma_cP
    [Sigma_cP,logPost] = imhSigma_cP(Sigma_cP,...
                                     parsSigma_cP_mle,...
                                     covSigma_cP_mle,...
                                     EX_Sigma_cP,...
                                     VAR_Sigma_cP,...
                                     EX_Qonly,...
                                     VAR_Qonly,...
                                     Y,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,sige,...
                                     fi);
    
    % ### independence Metropolis-Hastings step for kinfQ and K1Q_X
    [kinfQ,K1Q_X,logPost] = imhQonly(logPost,...
                                     kinfQ,...
                                     K1Q_X,...
                                     parsQonly_mle,...
                                     covQonly_mle,...
                                     EX_Qonly,...
                                     VAR_Qonly,...
                                     EX_Sigma_cP,...
                                     VAR_Sigma_cP,...
                                     Y,W,cP,mats,dt,lam0,lam1,Sigma_cP,sige,...
                                     fi);
                                 
    % P and Q persistence expressed via eigenvalues, implied yields calcs
    [~,K1Q_cP,AcP,BcP,AX,BX,AXL,BXL] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);
    I_K1Q_cP = K1Q_cP + IN;
    I_K1P_cP = I_K1Q_cP + lam1;
    
    out_eig_I_K1Q_cP(:,iter) = sort(eig(I_K1Q_cP),'descend');
    out_eig_I_K1P_cP(:,iter) = sort(eig(I_K1P_cP),'descend');
    out_AcP(:,iter)          = AcP;
    out_BcP(:,iter)          = reshape(BcP,J*N,1);
    out_AX(:,iter)           = AX;
    out_BX(:,iter)           = reshape(BX,J*N,1);
    out_AXL(:,iter)          = AXL;
    out_BXL(:,iter)          = reshape(BXL,JL*N,1);
    
    % store posterior draws
    if isempty(varargin)
        out_kinfQ(iter)      = kinfQ;                   
        out_K1Q_X(:,iter)    = diag(K1Q_X);            
        out_lam0(:,iter)     = lam0;
        out_lam1(:,iter)     = reshape(lam1,N*N,1);
        out_Sigma_cP(:,iter) = Sigma_cP2pars(Sigma_cP);
        out_sige(iter)       = sige;
        out_gamma(:,iter)    = gamma';
        
        out_gamma_prob(:,iter) = prob';
        
    elseif strcmp(varargin{1},'transform')
        out_kinfQ(iter)      = scalekinfQ(kinfQ);
        out_K1Q_X(:,iter)    = diag(mlgits(K1Q_X));
        out_lam0(:,iter)     = lam0;                % Gibbs step
        out_lam1(:,iter)     = reshape(lam1,N*N,1); % Gibbs step
        out_Sigma_cP(:,iter) = Sigma_cP2pars(logd(Sigma_cP));
        out_sige(iter)       = sige;                % Gibbs step
        out_gamma(:,iter)    = gamma';
        
        out_gamma_prob(:,iter) = prob';
        
    else
        out_kinfQ(iter)      = kinfQ;
        out_K1Q_X(:,iter)    = diag(K1Q_X);
        out_lam0(:,iter)     = lam0;
        out_lam1(:,iter)     = reshape(lam1,N*N,1);
        out_Sigma_cP(:,iter) = Sigma_cP2pars(Sigma_cP);
        out_sige(iter)       = sige;
        out_gamma(:,iter)    = gamma';
        
        out_gamma_prob(:,iter) = prob';
        
    end
    % store log posterior
    out_logPost(iter) = logPost;
                    
end

%% output for IBIS
out = [out_logPost',out_kinfQ',out_K1Q_X',out_lam0',out_lam1',out_Sigma_cP',out_sige'];

%% acceptance rates
acr  = acrate([out_Sigma_cP(1,:);out_kinfQ]);

acr_ = struct();

acr_.acr_Sigma_cP    = acr(1);
acr_.acr_kinfQ_K1Q_X = acr(2);

%% full output
out_ = struct();

out_.out_eig_I_K1Q_cP = out_eig_I_K1Q_cP;
out_.out_eig_I_K1P_cP = out_eig_I_K1P_cP;
out_.out_AcP          = out_AcP;
out_.out_BcP          = out_BcP;
out_.out_AX           = out_AX;
out_.out_BX           = out_BX;
out_.out_AXL          = out_AXL;
out_.out_BXL          = out_BXL;

out_.out_kinfQ    = out_kinfQ;
out_.out_K1Q_X    = out_K1Q_X;
out_.out_lam0     = out_lam0;
out_.out_lam1     = out_lam1;
out_.out_Sigma_cP = out_Sigma_cP;
out_.out_sige     = out_sige;

out_.out_gamma    = out_gamma;
out_.out_gamma_prob = out_gamma_prob;

out_.out_logPost  = out_logPost;

end