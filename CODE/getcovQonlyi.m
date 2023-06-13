function covQonly = getcovQonlyi(eps,lift,kinfQ,K1Q_X,Y,W,cP,mats,dt,lam0,lam1,Sigma_cP,sige)
% This function calculates covariance matrix for Q parameters at the ML 
% estimate as inverse of negative hessian (Fisher information matrix)
% and increments its diagonal elements if specified

% stack Q parameters into a vector, scale and transform
parsQonly = Qonly2pars(scalekinfQ(kinfQ),mlgits(K1Q_X));

% get hessian at Q parameters
hessQonly = numhessi(@hessLlkQonly,parsQonly,eps,Y,W,cP,mats,dt,lam0,lam1,Sigma_cP,sige);

% calculate covariance matrix as minus inverse of the Hessian
covQonly = makepd(hessQonly\eye(numel(parsQonly)));

% increment diagonal elements
idx = (eye(numel(parsQonly)) == 1);
covQonly(idx) = (1 + lift)*covQonly(idx);
covQonly = makepd(covQonly);

end