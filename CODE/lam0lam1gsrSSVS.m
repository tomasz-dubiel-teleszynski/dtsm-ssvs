function [lam0,lam1] = lam0lam1gsrSSVS(lam0,lam1,...
                                       gamma,...
                                       EX_lam0lam1,...
                                       VAR_lam0lam1,...
                                       W,cP,mats,dt,kinfQ,K1Q_X,Sigma_cP,...
                                       varargin)
% This function implements Gibbs sampler for lam0 and lam1 with
% non-zero restrictions assuming joint MVN(lam_,V_) prior

% OLS/MLE/GLS for lam0 and lam1
% get Q parameters for observed factors cP from observed factor loading
% matrix and Q parameters for unobserved factors X
[K0Q_cP,K1Q_cP] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);

% get number of observations and number of factors
[T,N] = size(cP);

% tempering
if isempty(varargin)
    fi = 1;
else
    fi = varargin{1};
end

% vectorize
X = cP';
vecX = reshape(X(:,2:T),N*(T-1),1);
Z = [ones(1,T-1); X(:,1:T-1)];

% tempering step -------------------------------
vecX(end-N+1:end) = sqrt(fi)*vecX(end-N+1:end);
Z(:,end) = sqrt(fi)*Z(:,end);
% ----------------------------------------------

gamma0 = ones(1,N*(N + 1));

INN1 = eye(N*(N + 1));
idx = (gamma0 == 1);
S = INN1(:,idx);
r = [K0Q_cP;reshape(K1Q_cP + eye(N),N*N,1)];
IN = eye(N);
z = vecX - kron(Z',IN)*r;

% GLS estimate
Omega_cP = Sigma_cP*Sigma_cP';
iOmega_cP = Omega_cP\IN;
nzr = sum(gamma0);
INZR = eye(nzr);
Vhat = (S'*kron(Z*Z',iOmega_cP)*S)\INZR;
lamhat = Vhat*(S'*kron(Z,iOmega_cP)*z);

% SSVS prior
tau = VAR_lam0lam1.tau0.*(1 - gamma') + VAR_lam0lam1.tau1.*gamma';
iDRD = diag(1./(tau.*tau)); 

% Gibbs sampler
Vbar = (iDRD + Vhat\INZR)\INZR;
lambar = Vbar*(Vhat\lamhat);
% Cholesky decomposition
L = chol(Vbar,'lower');
lam0lam1 = lambar + L*randn(nzr,1);

% check if transition matrix under P for observed factors is not
% explosive that is check if its eigenvalues are below 1 in absolute
% value
[~,K1Q_cP] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);
tmp  = gamma0;
tmp(idx) = lam0lam1;
tmp1 = eig(K1Q_cP + reshape(tmp(N+1:end),N,N) + IN);
n = 1;
while any(abs(tmp1) >= 1) && (n <= 1)
    lam0lam1 = lambar + L*randn(nzr,1);
    tmp  = gamma0;
    tmp(idx) = lam0lam1;
    tmp1 = eig(K1Q_cP + reshape(tmp(N+1:end),N,N) + IN);
    n = n + 1;
end

% proposed lam0 and lam1
if all(abs(tmp1) < 1)
    tmp  = gamma0;
    tmp(idx) = lam0lam1;   
    lam0 = tmp(1:N)';
    lam1 = reshape(tmp(N+1:end),N,N);
end

end