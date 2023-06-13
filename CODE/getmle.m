function [kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige] = getmle(gamma,...
                                                        Y0,W,cP0,mats,dt,...
                                                        K1Q_X0,...
                                                        Sigma_cP0,...
                                                        satf)

% structure with additional variables
rest = struct();
rest.gamma = gamma;
rest.Y = Y0;
rest.W = W;
rest.cP = cP0;
rest.mats = mats;
rest.dt = dt;                                                                            

% starting values
pars0 = theta2pars(mlgitsmle(K1Q_X0),logdmle(Sigma_cP0));

% optimization 
pars = getoptim(pars0,rest,satf);

% parameter estimates
[K1Q_X,Sigma_cP] = pars2theta(pars);
K1Q_X = mlgitinvsmle(K1Q_X);
Sigma_cP = expdmle(Sigma_cP);

% remaining parameter estimates for kinfQ, sige, lam0 and lam1
[Llk,kinfQ,sige,lam0,lam1] = jszllk(gamma,Y0,W,cP0,mats,dt,[],K1Q_X,Sigma_cP,[]);

Llk,kinfQ,K1Q_X,lam0,lam1,Sigma_cP,sige

end