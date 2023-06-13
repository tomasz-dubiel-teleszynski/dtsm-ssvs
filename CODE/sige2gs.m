function sige2 = sige2gs(alpha_sige2,...
                           beta_sige2,...
                           Y,W,cP,mats,dt,kinfQ,K1Q_X,Sigma_cP,...
                           varargin)
% This function implements Gibbs sampler for sige2 assuming IG(alpha,beta)
% prior

% get Ricatti loadings for observed yields from factor loading matrix and
% Q parameters for unobserved factors X
[~,~,AcP,BcP] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);

% get number of observations
T = size(Y,1);

% get number of factors and number of yields
[N,J] = size(W);

% tempering
if isempty(varargin)
    fi = 1;
else
    fi = varargin{1};
end

% Gibbs sampler
Y_hat = repmat(AcP',T,1) + cP*BcP';
%errors = (Y - Y_hat)*null(W); % null space of W
errors = Y - Y_hat;

% tempering step -----------------------
errors(end,:) = sqrt(fi)*errors(end,:);
% --------------------------------------

% proposed sige2
sige2 = 1./gamrnd((alpha_sige2 + (T - 1 + fi)*(J - N))/2,1/((beta_sige2 + sum(diag(errors*errors')))/2));

end