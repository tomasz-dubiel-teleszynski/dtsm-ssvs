function [kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige] = pars2svms(pars,dims)
% This function distributes all parameters stacked into vector pars into
% scalars, vectors and matrices in the model 

% structure for dimensions of all parameters from vectors and matrices
% in the model 
d = dims;

% parameter vector to scalars, vectors and matrices
% kinfQ
idx1 = 0; idx2 = idx1 +1;
kinfQ = reshape(pars(idx1+1:idx2),1,1);

% K1Q_X
N = d.nrlam0;
tmp = d.nrK1Q_X*d.ncK1Q_X-(N*N - N);
idx1 = idx2; idx2 = idx2 +tmp;
K1Q_X = diag(pars(idx1+1:idx2));

% lam0
idx1 = idx2; idx2 = idx2 +d.nrlam0*d.nclam0;
lam0 = reshape(pars(idx1+1:idx2),d.nrlam0,d.nclam0);

% lam1
idx1 = idx2; idx2 = idx2 +d.nrlam1*d.nclam1;
lam1 = reshape(pars(idx1+1:idx2),d.nrlam1,d.nclam1);

% Sigma_cP
tmp = d.nrSigma_cP*d.ncSigma_cP-(N*N - N)/2;
idx1 = idx2; idx2 = idx2 +tmp;
Sigma_cP = pars(idx1+1:idx2);
if N == 3
    Sigma_cP = [Sigma_cP(1:d.nrSigma_cP);0;Sigma_cP(d.nrSigma_cP+1:tmp-1);0;0;Sigma_cP(tmp)];
elseif N == 2
    Sigma_cP = [Sigma_cP(1:d.nrSigma_cP);0;Sigma_cP(tmp)];
end
Sigma_cP = reshape(Sigma_cP,d.nrSigma_cP,d.ncSigma_cP);

% sige
idx1 = idx2; idx2 = idx2 +1;
sige = pars(idx1+1:idx2);

end