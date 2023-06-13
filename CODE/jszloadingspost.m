function [K0Q_cP,K1Q_cP,AcP,BcP,AX,BX,AXL,BXL,AcPL,BcPL] = jszloadingspost(W,mats,dt,kinfQ,K1Q_X,Sigma_cP)
% This function calculates Q parameters for observed factors cP and Ricatti
% loading matrices for observed yields from real factor loading matrix and
% Q parameters for unobserved factors X

% get number of factors
N = size(K1Q_X,1);

% scaled maturities
mats_periods = round(mats/dt);

% adjust K1Q_X
[K1Q_X,~,m1] = jszAdjustK1QX(K1Q_X);

%  Q drift in VAR(1) for unobserved factors X
K0Q_X     = zeros(N,1);
K0Q_X(m1) = kinfQ;

% deterministic terms
rho0d = 0;
rho1d = ones(N,1);

% preliminary iteration of Ricatti equations
[~,BX,~,BXL] = gaussianloadingspost(mats_periods,dt,K0Q_X,K1Q_X,zeros(N),rho0d*dt,rho1d*dt);

% extract covariance matrix of VAR(1) for unobserved factors from
% covariance matrix of VAR(1) for observed factors 
WBX = W*BX;
Omega_cP = Sigma_cP*Sigma_cP';
Omega_X  = (WBX\Omega_cP)/WBX'; % equivalent to: Omega_X = inv(WBX)*Omega_cP*inv(WBX)';

% final iteration of Ricatti equations
[AX,~,AXL] = gaussianloadingspost(mats_periods,dt,K0Q_X,K1Q_X,Omega_X,rho0d*dt,rho1d*dt);

% Q parameters for observed factors cP 
WAX = W*AX;
K1Q_cP = WBX*(K1Q_X/WBX); % equivalent to: K1Q_cP = WBX*K1Q_X*inv(WBX);
K0Q_cP = WBX*K0Q_X - K1Q_cP*WAX;

% rotated Ricatti loadings 
BcP = BX/WBX;          % equivalent to: BcP = BX*inv(WBX);
AcP = AX - BcP*WAX;    % equivalent to: AcP = AX - BX*inv(WBX)*WAX; 

% rotated large Ricatti loadings
BcPL = BXL/WBX;        % equivalent to: BcPL = BXL*inv(WBX);
AcPL = AXL - BcPL*WAX; % equivalent to: AcPL = AXL - BXL*inv(WBX)*WAX;

end