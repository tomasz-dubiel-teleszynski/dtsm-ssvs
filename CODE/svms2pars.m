function [pars,dims] = svms2pars(kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige)
% This function stacks all parameters from scalars, vectors and matrices
% in the model into a column vector

% structure for dimensions of all parameters from vectors and matrices
% in the model 
d = struct();

% K1Q_X
[d.nrK1Q_X,d.ncK1Q_X] = size(K1Q_X); 
K1Q_X = diag(K1Q_X);

% lam0
[d.nrlam0,d.nclam0] = size(lam0);

% lam1
[d.nrlam1,d.nclam1] = size(lam1);

% Sigma_cP
[d.nrSigma_cP,d.ncSigma_cP] = size(Sigma_cP);
N = numel(lam0);
if N == 3
    Sigma_cP = [Sigma_cP(:,1);Sigma_cP(2:3,2);Sigma_cP(3,3)];
elseif N == 2
     Sigma_cP = [Sigma_cP(:,1);Sigma_cP(2,2)];
end
    
% scalars, vectors and matrices stacked into a vector
pars =[kinfQ;...
       reshape(K1Q_X,d.nrK1Q_X*d.ncK1Q_X-(N*N - N),1);... 
       reshape(lam0,d.nrlam0*d.nclam0,1);...
       reshape(lam1,d.nrlam1*d.nclam1,1);...
       reshape(Sigma_cP,d.nrSigma_cP*d.ncSigma_cP-((N*N - N)/2),1);... 
       sige];
dims = d;
   
end