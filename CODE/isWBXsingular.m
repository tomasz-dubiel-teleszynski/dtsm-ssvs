function [issingular, WBX] = isWBXsingular(W,mats,dt,kinfQ,K1Q_X)
% This function checks if product of matrices W and BX where the latter
% stems from Ricatti equations is singular or not

% get number of factors
N = size(K1Q_X,1);

% scaled maturities
mats_periods = round(mats/dt);

%  Q drift in VAR(1) for unobserved factors X
K0Q_X = zeros(N,1);
K0Q_X(1,1) = kinfQ;

% deterministic terms
rho0d = 0;
rho1d = ones(N,1);

% preliminary iteration of Ricatti equations
[~,BX] = gaussianloadings(mats_periods,dt,K0Q_X,K1Q_X,zeros(N),rho0d*dt,rho1d*dt);

% extract covariance matrix of VAR(1) for unobserved factors from
% covariance matrix of VAR(1) for observed factors 
WBX = W*BX;

% check singularity of WBX
c = cond(WBX);
if c > 1e15
    issingular = true;
else
    issingular = false;
end