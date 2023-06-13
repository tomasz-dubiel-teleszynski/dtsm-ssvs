function lamsd = gprior(g,Y0,W,cP0,mats,dt)

% get number of observations and number of factors
[T,N] = size(cP0); 

% get non-zeros restrictions
model = 'MAXIMALLY_FLEXIBLE';
gamma = gammarest(model,N);

% simulated annealing true of false
satf = false;

% number of seeds
nSeeds = 10000;

% OLS/ML estimate for Sigma_cP from observed factors cP
Sigma_cP0 = getolsSigma_cP(cP0);

% get starting values for K1Q_X
K1Q_X0 = getstartingvaluesformle(nSeeds,gamma,Y0,W,cP0,mats,dt,Sigma_cP0);

% get MLEs
[~,~,~,~,Sigma_cP_mle] = getmle(gamma,...
                                Y0,W,cP0,mats,dt,...
                                K1Q_X0,...
                                Sigma_cP0,...
                                satf);
                     
% vectorize
X = cP0';
Z = [ones(1,T-1); X(:,1:T-1)];

INN1 = eye(N*(N + 1));
idx = (gamma == 1);
S = INN1(:,idx);
IN = eye(N);

% GLS estimate
Omega_cP = Sigma_cP_mle*Sigma_cP_mle';
iOmega_cP = Omega_cP\IN;
nzr = sum(gamma);
INZR = eye(nzr);
Vhat = (S'*kron(Z*Z',iOmega_cP)*S)\INZR;

% standard deviations for lam0 and lam1 in diagonalized g-prior
lamsd = sqrt(g)*sqrt(diag(Vhat));

end

