function [K0Q_cP, K1Q_cP, rho0_cP, rho1_cP] = jszrotation(W, K1Q_X, K0Q_X, rho0_X, rho1_X, BX, AX)

WBX = W*BX;
WAX = W*AX;

K1Q_cP = WBX*(K1Q_X/WBX); % equivalent to: K1Q_cP = WBX*K1Q_X*inv(WBX);
K0Q_cP = WBX*K0Q_X - K1Q_cP*WAX;

rho0_cP = rho0_X - rho1_X'*(WBX\WAX);
rho1_cP = (WBX')\rho1_X;
