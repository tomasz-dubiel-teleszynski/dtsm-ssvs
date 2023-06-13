function Llk = objmle(pars,rest) 

gamma = rest.gamma;
Y = rest.Y;
W = rest.W;
cP = rest.cP;
mats = rest.mats;
dt = rest.dt;

[K1Q_X,Sigma_cP]= pars2theta(pars);
K1Q_X = mlgitinvsmle(K1Q_X);
Sigma_cP = expdmle(Sigma_cP);

Llk = jszllk(gamma,Y,W,cP,mats,dt,[],K1Q_X,Sigma_cP,[]);