function [Llk,kinfQ,sige,lam0,lam1] = jszllk(gamma,Y,W,cP,mats,dt,kinfQ,K1Q_X,Sigma_cP,sige)

% negative Q likelihood

J = numel(mats);
[T,N] = size(cP(2:end,:));

if isempty(kinfQ)
    rho0_cP = 0;
    [K1Q_X,~,~,~,BcP,~,BX,alpha0_cP,alpha1_cP,alpha0_X,alpha1_X,~,~,m1] = jszloadingsrho0cP(W,mats,dt,rho0_cP,K1Q_X,Sigma_cP);
    V = null(W)';
    kinfQ = ((mean(Y(2:end,:))' - alpha1_cP - BcP*mean(cP(2:end,:)).')'*((V')*V*alpha0_cP))/(alpha0_cP'*(V')*V*alpha0_cP);
    AX = alpha0_X*kinfQ + alpha1_X;
    AcP = alpha0_cP*kinfQ + alpha1_cP;
    K0Q_X = zeros(N,1);
    K0Q_X(m1) = kinfQ;
    [K0Q_cP, K1Q_cP] = jszrotation(W, K1Q_X, K0Q_X, 0, ones(N,1), BX, AX);
else
    [K0Q_cP,K1Q_cP,AcP,BcP] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);
end
Y_hat = repmat(AcP',T,1) + cP(2:end,:)*BcP'; 
errors = Y(2:end,:) - Y_hat;

if isempty(sige)
    sige = sqrt(sum(errors(:).^2)/(T*(J - N)));
end
    
llkQ = - (J - N)*0.5*log(2*pi) - (J - N)*log(sige) - 0.5*diag(errors*errors')/(sige*sige);
LlkQ = sum(llkQ);
LlkQ = - LlkQ;

% negative P likelihood

T = size(cP,1);
X = cP';
vecX = reshape(X(:,2:T),N*(T-1),1);
Z = [ones(1,T-1); X(:,1:T-1)];
INN1 = eye(N*(N + 1));
idx = (gamma == 1);
S = INN1(:,idx);
r = [K0Q_cP;reshape(K1Q_cP + eye(N),N*N,1)];
IN = eye(N);
z = vecX - kron(Z',IN)*r;

Omega_cP = Sigma_cP*Sigma_cP';
iOmega_cP = Omega_cP\IN;
nzr = sum(gamma);
Vhat = (S'*kron(Z*Z',iOmega_cP)*S)\eye(nzr);
lamhat = Vhat*(S'*kron(Z,iOmega_cP)*z);

tmp  = gamma;
tmp(idx) = lamhat;

lam0 = tmp(1:N)';
lam1 = reshape(tmp(N+1:end),N,N);

[K0P_cP,I_K1P_cP] = getPdyn(K0Q_cP,K1Q_cP,lam0,lam1);
innovations = cP(2:T,:) - (repmat(K0P_cP',T-1,1) + cP(1:T-1,:)*I_K1P_cP');

llkP = - 0.5*N*log(2*pi) - 0.5*log(det(Omega_cP)) - 0.5*diag((innovations/Omega_cP)*innovations'); 
LlkP = sum(llkP);
LlkP = - LlkP;

% negative likelihood

Llk = LlkQ + LlkP;

end

