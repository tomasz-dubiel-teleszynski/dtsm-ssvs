function [A,B,AL,BL] = gaussianloadings(mats_periods,dt,K0d,K1d,Omega_d,rho0d,rho1d)
% This function calculates Ricatti loadings

% get number of factors
N = length(K0d);
% get number of yields
J = length(mats_periods);
% to annualize maturities
n_per = 12;

% initial values for Ricatti loadings
Atemp = 0;
Btemp = zeros(N,1);

% set containers for Ricatti loadings
Ay = zeros(1,J);
By = zeros(N,J);

% set containers for large Ricatti loadings
mats_periodsL = 12:n_per:mats_periods(end);
JL = numel(mats_periodsL);

AyL = zeros(1,JL);
ByL = zeros(N,JL);

% transpose
K0dp = K0d';
K1dp = K1d';

% iterate Ricatti equations
curr_mat = 1;
curr_matL = 1;
for i = 1:mats_periods(J)

    % Ricatti equations
    Atemp = Atemp + K0dp*Btemp + 0.5*Btemp'*Omega_d*Btemp - rho0d;
    Btemp = Btemp + K1dp*Btemp - rho1d;
    
    % select only chosen maturities
    if i == mats_periods(curr_mat)
        Ay(1,curr_mat) = -Atemp/mats_periods(curr_mat);
        By(:,curr_mat) = -Btemp/mats_periods(curr_mat);
        curr_mat = curr_mat + 1;
    end
    
    % select all yearly maturities
    if i == mats_periodsL(curr_matL)
        AyL(1,curr_matL) = -Atemp/mats_periodsL(curr_matL);
        ByL(:,curr_matL) = -Btemp/mats_periodsL(curr_matL);
        curr_matL = curr_matL + 1;
    end
end

% scale Ricatti loadings with time step
A = Ay'/dt;
B = By'/dt;

% scale large Ricatti loadings with time step
AL = AyL'/dt;
BL = ByL'/dt;

end