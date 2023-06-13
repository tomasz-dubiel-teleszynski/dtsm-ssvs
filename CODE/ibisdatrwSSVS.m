function ibisdatrwSSVS()

%% timer on
tic;

%% fix seed
seed = 616;
rng(seed);

%% data specification
% set number of factors
N = 3;
% monthly maturities
mats = [12 24 36 48 60 84 120]';
% time step one month
dt = 1;

%% models
model = 'MAXIMALLY_FLEXIBLE';

% for MLE and MCMC to run
exists_MLE = false; exists_MCMC = false;        

%% restrictions on lam0 and lam1
% get non-zeros restrictions
gamma = gammarest(model,N);
% number of non-zero restrictions
nzr = sum(gamma);

%% use real yields and extract observed factors cP from real yields
start_date = 19900101; %19850101
end_date   = 20190000; %20080000
cut_date   = 20080000; %19970000
fct_date   = 20080000; %19970000
% related strings
start_year_str = num2str(start_date);    start_year_str = start_year_str(1:4);
end_year_str   = num2str(end_date - 1e4);  end_year_str = end_year_str(1:4);    
cut_year_str   = num2str(cut_date - 1e4);  cut_year_str = cut_year_str(1:4);
fct_year_str   = num2str(fct_date - 1e4);  fct_year_str = fct_year_str(1:4);

% real yields and observed factors at maturities in mats
% true out-of-sample
[~,Y,~,dates] = realWYcP(mats,N,start_date,end_date);
W = realWYcP(mats,N,start_date,fct_date);
cP = Y*W';

% real yields at maturities from 1 to mats(end)
matsL = 1:mats(end);
[~,YL] = realWYcP(matsL,N,start_date,end_date);

% shift variable
shift_to_obs = 12;

% get number of observations
T = size(Y(1:end-shift_to_obs,:),1);

%% cut yields and factors to first half of the observations 
% this is done in order to initiate prior sample of particles
cut = find(dates < cut_date,1,'last'); 
fct = find(dates < fct_date,1,'last'); 

% for MCMC
Y0  = Y(1:cut,:);
cP0 = cP(1:cut,:);

% for g-prior
Yg  = Y(1:fct,:);
cPg = cP(1:fct,:);

%% priors
% univariate IG(alpha,beta) prior for sige2
alpha_sige2 = 0;
beta_sige2  = 0;

% diagonalized g-prior for lam0 and lam1
g = T;
lamsd = gprior(g,Yg,W,cPg,mats,dt);
VAR_lam0lam1 = diag(lamsd(logical(gamma)).^2);

% SSVS prior
c1 = sqrt(g);
c0 = 1/c1;
tau0 = c0*sqrt(diag(VAR_lam0lam1/g));
tau1 = c1*tau0/c0;
EX_lam0lam1 = 0.5*ones(nzr,1); % Binomial(numel(gamma),0.5)
VAR_lam0lam1 = struct();
VAR_lam0lam1.tau0 = tau0;
VAR_lam0lam1.tau1 = tau1;

% univariate N(EX_,VAR_) priors for each element of Sigma_cP
EX_Sigma_cP  = 0;
VAR_Sigma_cP = 1e6;

% univariate N(EX_,VAR_) priors for kinfQ and each element of K1Q_X
EX_Qonly  = 0;
VAR_Qonly = 1e6;

%% initial MLEs

if ~exists_MLE
    % simulated annealing true of false
    satf = false;

    % number of seeds
    nSeeds = 10000;

    % OLS/ML estimate for Sigma_cP from observed factors cP
    Sigma_cP0 = getolsSigma_cP(cP0);

    % get starting values for K1Q_X
    K1Q_X0 = getstartingvaluesformle(nSeeds,gamma,Y0,W,cP0,mats,dt,Sigma_cP0);

    % get MLEs
    [kinfQ_mle,K1Q_X_mle,lam0_mle,lam1_mle,Sigma_cP_mle,sige_mle] = getmle(gamma,...
                                                                           Y0,W,cP0,mats,dt,...
                                                                           K1Q_X0,...
                                                                           Sigma_cP0,...
                                                                           satf);
end

% one can also comment above MLE part and just load MLEs from a run made 
% in advance
if exists_MLE
    eval(['load mles_ibisd_',fct_year_str,'_SSVS;']);
else
    eval(['save mles_ibisd_',fct_year_str,'_SSVS kinfQ_mle K1Q_X_mle lam0_mle lam1_mle Sigma_cP_mle sige_mle;']);
end   

% transform and scale Sigma_cP_mle for independence Metropolis-Hastings sampler
parsSigma_cP_mle = Sigma_cP2pars(logd(Sigma_cP_mle));

% transform and scale kinfQ_mle and K1Q_X_mle for independence Metropolis-Hastings sampler
parsQonly_mle = Qonly2pars(scalekinfQ(kinfQ_mle),mlgits(K1Q_X_mle));

%% parameter covariance matrices at initial MLEs
[covSigma_cP_mle,covQonly_mle] = getcovsi(kinfQ_mle,...
                                          K1Q_X_mle,...
                                          lam0_mle,...
                                          lam1_mle,...
                                          Sigma_cP_mle,...
                                          sige_mle,...
                                          Y0,W,cP0,mats,dt);

%% MCMC output for IBIS
% this argument is normally provided within IBIS algorithm but here it is 
% initially filled with MLEs
[in,dims] = svms2pars(kinfQ_mle,...
                      K1Q_X_mle,...
                      lam0_mle,...
                      lam1_mle,...
                      Sigma_cP_mle,...
                      sige_mle);
% number of parameters
K = numel(in);
                            
%% inputs required in MCMC
% initial MCMC specification
Niter_init = 100000;
% take every n-th iteration and form initial sample of particles
n = 50;
idx = n:n:Niter_init;
% parallelize MCMC
Niter_init_mc = 11000;
burnin = 1000;
MC = 10;
OUT = cell(MC,1);
ACR_ = cell(MC,1);
OUT_ = cell(MC,1);

%% initialization

if ~exists_MCMC
    % run initial MCMC using first half of the observations to get initial
    % sample of particles
    parfor mc = 1:MC
        [OUT{mc},ACR_{mc},OUT_{mc}] = mhwg2ibisiSSVS(Niter_init_mc,...
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
                                                     Y0,W,cP0,mats,dt,...
                                                     'transform');        
    end
    
    for mc = 1:MC
        if mc > 1
            out = [out;OUT{mc}(burnin+1:end,:)];  
            acr_.acr_Sigma_cP = [acr_.acr_Sigma_cP,ACR_{mc}.acr_Sigma_cP];  
            acr_.acr_kinfQ_K1Q_X = [acr_.acr_kinfQ_K1Q_X,ACR_{mc}.acr_kinfQ_K1Q_X];
            
            out_.out_eig_I_K1Q_cP = [out_.out_eig_I_K1Q_cP,OUT_{mc}.out_eig_I_K1Q_cP(:,burnin+1:end)];
            out_.out_eig_I_K1P_cP = [out_.out_eig_I_K1P_cP,OUT_{mc}.out_eig_I_K1P_cP(:,burnin+1:end)];
            out_.out_AcP          = [out_.out_AcP,OUT_{mc}.out_AcP(:,burnin+1:end)];
            out_.out_BcP          = [out_.out_BcP,OUT_{mc}.out_BcP(:,burnin+1:end)];
            out_.out_AX           = [out_.out_AX,OUT_{mc}.out_AX(:,burnin+1:end)];
            out_.out_BX           = [out_.out_BX,OUT_{mc}.out_BX(:,burnin+1:end)];
            out_.out_AXL          = [out_.out_AXL,OUT_{mc}.out_AXL(:,burnin+1:end)];
            out_.out_BXL          = [out_.out_BXL,OUT_{mc}.out_BXL(:,burnin+1:end)];
            out_.out_kinfQ        = [out_.out_kinfQ,OUT_{mc}.out_kinfQ(:,burnin+1:end)];
            out_.out_K1Q_X        = [out_.out_K1Q_X,OUT_{mc}.out_K1Q_X(:,burnin+1:end)];
            out_.out_lam0         = [out_.out_lam0,OUT_{mc}.out_lam0(:,burnin+1:end)];
            out_.out_lam1         = [out_.out_lam1,OUT_{mc}.out_lam1(:,burnin+1:end)];
            out_.out_Sigma_cP     = [out_.out_Sigma_cP,OUT_{mc}.out_Sigma_cP(:,burnin+1:end)];
            out_.out_sige         = [out_.out_sige,OUT_{mc}.out_sige(:,burnin+1:end)];
            
            out_.out_gamma        = [out_.out_gamma,OUT_{mc}.out_gamma(:,burnin+1:end)];
            out_.out_gamma_prob   = [out_.out_gamma_prob,OUT_{mc}.out_gamma_prob(:,burnin+1:end)];
            
            out_.out_logPost      = [out_.out_logPost,OUT_{mc}.out_logPost(:,burnin+1:end)]; 
            
        else
            out = OUT{mc}(burnin+1:end,:);
            acr_ = ACR_{mc};
            
            out_ = OUT_{mc};
            out_.out_eig_I_K1Q_cP = out_.out_eig_I_K1Q_cP(:,burnin+1:end);
            out_.out_eig_I_K1P_cP = out_.out_eig_I_K1P_cP(:,burnin+1:end);
            out_.out_AcP          = out_.out_AcP(:,burnin+1:end);
            out_.out_BcP          = out_.out_BcP(:,burnin+1:end);
            out_.out_AX           = out_.out_AX(:,burnin+1:end);
            out_.out_BX           = out_.out_BX(:,burnin+1:end);
            out_.out_AXL          = out_.out_AXL(:,burnin+1:end);
            out_.out_BXL          = out_.out_BXL(:,burnin+1:end);
            out_.out_kinfQ        = out_.out_kinfQ(:,burnin+1:end);
            out_.out_K1Q_X        = out_.out_K1Q_X(:,burnin+1:end);
            out_.out_lam0         = out_.out_lam0(:,burnin+1:end);
            out_.out_lam1         = out_.out_lam1(:,burnin+1:end);
            out_.out_Sigma_cP     = out_.out_Sigma_cP(:,burnin+1:end);
            out_.out_sige         = out_.out_sige(:,burnin+1:end);
            
            out_.out_gamma        = out_.out_gamma(:,burnin+1:end);
            out_.out_gamma_prob   = out_.out_gamma_prob(:,burnin+1:end);
            
            out_.out_logPost      = out_.out_logPost(:,burnin+1:end);
            
        end 
    end
    
end

% one can also comment out this function above and just load particles from
% MCMC run of Niter_init iterations made in advance
if exists_MCMC
    eval(['load particles_ibisd_',fct_year_str,'_',num2str(Niter_init),'_',num2str(n),'_SSVS;']);
else
    eval(['save particles_ibisd_',fct_year_str,'_',num2str(Niter_init),'_',num2str(n),'_SSVS out acr_ out_ idx;']);
end   
                                
% select every 50th sample as initial particles
particles       = out(idx,2:end);

indicators      = out_.out_gamma(:,idx);
indicators_prob = out_.out_gamma_prob(:,idx);

%% main algorithm
% particle parameters
Np   = size(particles,1); % number of particles
gam  = 0.7;               % percentage of Np particles to trigger resampling/rejuvenation
crit = gam*Np;            % ESS degeneracy criterion

% MCMC iterations in rejuvenation step
Niter = 5;

% set containers
us     =  cell(1,T - cut);
logus  =  cell(1,T - cut);
Ls     = zeros(1,T - cut);
moves  = zeros(1,T - cut);
jcs    = zeros(K,T - cut);
AVGACR = zeros(2,T - cut);

ws     =  cell(1,T - cut + 1);
logws  =  cell(1,T - cut + 1);
thetas =  cell(1,T - cut + 1);
ESSs   = zeros(1,T - cut + 1) + Np;

gammas      = cell(1,T - cut + 1);
gammas_prob = cell(1,T - cut + 1);

% initialization 
thetas{1} = particles;
ws{1}     = ones(Np,1);
logws{1}  = zeros(Np,1);

gammas{1} = indicators;
gammas_prob{1} = indicators_prob;

% containers
RXF = zeros(T - fct + 1,72);
RXO = RXF;
RXOBAR = RXF;
RXFALL = cell(T - fct + 1,1);
RXOALL = RXFALL;

% IBIS loop
for t = cut+1:T+1
    
    if t > fct
        [RXF(t - fct,:),RXO(t - fct,:),RXOBAR(t - fct,:),RXL,RXFALL{t - fct},RXOALL{t - fct}] = ...
            predablty(t - 1,ws{t - cut},thetas{t - cut},dims,YL,matsL,W,cP,mats,dt);
    end
    
    % terminate for loop after last 
    if t == T+1
       break;
    end
    
    % cut yields and factors to first t observations
    Yt  =  Y(1:t,:);
    cPt = cP(1:t,:);
    
    % truncate yields and factors to observation t and from t-1 to t respectively 
    Y1  =  Y(t,:);
    cP2 = cP(t-1:t,:);
    
    % select parameter container at time t
    thetast = thetas{t - cut};
    
    % incremental weights
    logiw = zeros(Np,1);
    parfor ip = 1:Np
        
        [kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige] = pars2svms(thetast(ip,:)',dims);
        kinfQ = rescalekinfQ(kinfQ);
        K1Q_X = mlgitinvs(K1Q_X);
        Sigma_cP = expd(Sigma_cP);
        
        logiw(ip) = getLlk1(Y1,W,cP2,mats,dt,...
                            kinfQ,...
                            K1Q_X,...
                            lam0,...
                            lam1,...
                            Sigma_cP,...
                            sige);
                        
    end
    
    % numerical fix
    logiw = numfix(logiw);
    
    % incremental weights
    logus{t - cut} = logiw;
    us{t - cut}    = exp(logus{t - cut});
    
    % update impotance weights
    logws{t - cut + 1} = logws{t - cut} + logus{t - cut};
    ws{t - cut + 1}    = exp(logws{t - cut +1});
    
    % effective sample size
    ESS = getESSnr(1,0,logws{t - cut + 1});
    % NaN and Inf checks
    if isnan(ESS) || isinf(ESS)
        ESS = 1;
    end
    
    if ESS > crit
        
        % and their weighted average
        Ls(t - cut) = sum(ws{t - cut}.*us{t - cut})/sum(ws{t - cut});
        
        thetas{t - cut + 1} = thetas{t - cut};
        ESSs(t - cut + 1)   = ESS;
        
        gammas{t - cut + 1}      = gammas{t - cut};
        gammas_prob{t - cut + 1} = gammas_prob{t - cut};
        
    else
        % record move
        moves(t - cut) = 1;
        
        % ------------------------------------------
        % START OF ADAPTIVE ITERMEDIATE TEMPERING --
        % ------------------------------------------
        
        % initalization
        r = 0;
        fi = 0;
        fi1 = 0;
        
        wnr = ws{t - cut};
        logwnr = logws{t - cut};
        Ls(t - cut) = 1;
        
        while fi < 1
           % counter 
           r = r + 1;
           
           if getESSnrat(1,fi1,logiw,logwnr) > crit
               % reset temperature
               fi = 1;
               
           else
               % compute temperature
               fi = bsESSnrat(fi1,logiw,logwnr,crit);
           
           end
           
           % incremental weights
           logus{t  - cut} = logiw*(fi - fi1);
           us{t - cut}    = exp(logus{t - cut});
           % and their weighted average
           Ls(t - cut) = Ls(t - cut)*sum(wnr.*us{t - cut})/sum(wnr);
           
           % compute unnormalized weights
           wnr = exp(logwnr + logiw*(fi - fi1));
           
           % BEFORE RESAMPLING
           if r > 1
               muhat    = sum(bsxfun(@times,wnr/sum(wnr),thetast1));
               tmp      = bsxfun(@minus,thetast1,muhat);
           else
               muhat    = sum(bsxfun(@times,wnr/sum(wnr),thetas{t - cut}));
               tmp      = bsxfun(@minus,thetas{t - cut},muhat);
           end
           Sigmahat = bsxfun(@times,wnr/sum(wnr),tmp)'*tmp;
           
           parsSigma_cP_mle = muhat(end-(N*N + N)/2:end-1)';
           parsQonly_mle = muhat(1:N+1)';
           covSigma_cP_mle = Sigmahat(end-(N*N + N)/2:end-1,end-(N*N + N)/2:end-1);
           covQonly_mle = Sigmahat(1:N+1,1:N+1);
           
           % RESAMPLING step using Svensson (2012)
           [~,j]  = histc(rand(Np,1), [0 cumsum((wnr/sum(wnr))')]);
           if r > 1
               thetar       = thetast1(j,:);
               gammasr      = gammast1(:,j);
           else
               thetar       = thetas{t - cut}(j,:);
               gammasr      = gammas{t - cut}(:,j);
           end
           
           % REJUVENATION step using MCMC kernel
           
           % select parameter container at time t+1
           thetast1 = thetas{t - cut + 1};
           gammast1      = gammas{t - cut + 1};
           gammas_probt1 = gammas_prob{t - cut + 1};
           
           % rejuvenation
           
           ACR = zeros(Np,2);
           parfor ip = 1:Np

               tmp = thetar(ip,:)';
               tmp(1) = rescalekinfQ(tmp(1));
               tmp(2:N+1) = diag(mlgitinvs(diag(tmp(2:N+1))));
               tmp(end-(N*N + N)/2:end-1) = Sigma_cP2pars(expd(pars2Sigma_cP(tmp(end-(N*N + N)/2:end-1))));

               % run MCMC using first t observations to rejuvenate particles
               [out,acr_,out_] = mhwg2ibisiSSVS(Niter,...
                                                tmp,... 
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
                                                gammasr(:,ip)',...
                                                Yt,W,cPt,mats,dt,...
                                                'transform',...
                                                fi);
               % store acceptance rates
               ACR(ip,:) = [acr_.acr_Sigma_cP,acr_.acr_kinfQ_K1Q_X];
               
               % replace old particle with rejuvenated one
               thetast1(ip,:)      = out(end,2:end);
               gammast1(:,ip)      = out_.out_gamma(:,end);
               gammas_probt1(:,ip) = out_.out_gamma_prob(:,end);
               
           end
           AVGACR(:,t - cut) = mean(ACR)';
           
           % reset unnormalized weights
           logwnr = zeros(Np,1);
           wnr = ones(Np,1);
           
           % intermediate incremental weights
           logiw = zeros(Np,1);
           parfor ip = 1:Np

               [kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige] = pars2svms(thetast1(ip,:)',dims);
               kinfQ = rescalekinfQ(kinfQ);
               K1Q_X = mlgitinvs(K1Q_X);
               Sigma_cP = expd(Sigma_cP);
               
               logiw(ip) = getLlk1(Y1,W,cP2,mats,dt,...
                                   kinfQ,...
                                   K1Q_X,...
                                   lam0,...
                                   lam1,...
                                   Sigma_cP,...
                                   sige);
                               
           end
           
           % numerical fix
           logiw = numfix(logiw);
           
           % replace temperatures
           fi1 = fi;
               
        end
        
        % ------------------------------------------
        % END OF ADAPTIVE ITERMEDIATE TEMPERING ----
        % ------------------------------------------
        
        % store after tempering
        thetas{t - cut + 1}      = thetast1;
        gammas{t - cut + 1}      = gammast1;
        gammas_prob{t - cut + 1} = gammas_probt1;
        
        % reset weights and ESS
        ws{t - cut + 1} = wnr;
        logws{t - cut + 1} = log(ws{t - cut + 1});
        ESSs(t - cut + 1) = getESSnr(1,fi1,logws{t - cut + 1});
        
        % correlations between non-jittered (thetar) and jittered (thetast1)
        % particles at each rejuvenation step
        jcsk = zeros(K,1);
        for k = 1:K
            jcsk(k) = corr(thetar(:,k),thetast1(:,k));
        end
        jcs(:,t - cut) = jcsk;
        
    end
    
    if t == T
        w        = ws{t - cut + 1}/sum(ws{t - cut +1});
        muhat    = sum(bsxfun(@times,w,thetas{t - cut + 1}));
        tmp      = bsxfun(@minus,thetas{t - cut + 1},muhat);
        Sigmahat = bsxfun(@times,w,tmp)'*tmp;
    end
    
end

% save results
eval(['save output_ibisdatr_',fct_year_str,'_SSVS;']);

%% timer off
toc

end