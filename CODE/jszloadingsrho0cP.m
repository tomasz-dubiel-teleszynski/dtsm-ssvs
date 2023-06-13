function [K1Q_X,K0Q_cP,K1Q_cP,AcP,BcP,AX,BX,alpha0_cP,alpha1_cP,alpha0_X,alpha1_X,rho0_cP,rho1_cP,m1] = jszloadingsrho0cP(W,mats,dt,rho0_cP,K1Q_X,Sigma_cP)

N = size(K1Q_X,1);
rho0d = 0;
rho1d = ones(N,1);
mats_periods = round(mats/dt);

[K1Q_X,~,m1] = jszAdjustK1QX(K1Q_X);

K0Q_X = zeros(N,1);
K0Q_X(m1) = 1;

[alpha0_X,BX] = gaussianloadings(mats_periods,dt,K0Q_X,K1Q_X,zeros(N),rho0d*dt,rho1d*dt);

WBX = W*BX;
Omega_cP = Sigma_cP*Sigma_cP';
Omega_X  = (WBX\Omega_cP)/WBX'; % equivalent to: Omega_X = inv(WBX)*Omega_cP*inv(WBX)';

AX1 = gaussianloadings(mats_periods,dt,K0Q_X,K1Q_X,Omega_X,rho0d*dt,rho1d*dt);

alpha1_X = AX1 - alpha0_X;

a0 = ones(1,N)*(WBX\W*alpha0_X);
a1 = ones(1,N)*(WBX\W*alpha1_X);

kinfQ = -(rho0_cP + a1)/a0;
K0Q_X(m1) = kinfQ;

AX = alpha0_X*kinfQ + alpha1_X;

alpha0_cP = alpha0_X - BX/WBX*(W*alpha0_X); 
alpha1_cP = alpha1_X - BX/WBX*(W*alpha1_X); 

WAX = W*AX;

BcP = BX/WBX;          % equivalent to: BcP = BX*inv(WBX);
AcP = AX - BcP*WAX;    % equivalent to: AcP = AX - BX*inv(WBX)*WAX; 

K1Q_cP = WBX*(K1Q_X/WBX); % equivalent to: K1Q_cP = WBX*K1Q_X*inv(WBX);
K0Q_cP = WBX*K0Q_X - K1Q_cP*WAX;

rho1_cP = WBX'\ones(N,1);

end