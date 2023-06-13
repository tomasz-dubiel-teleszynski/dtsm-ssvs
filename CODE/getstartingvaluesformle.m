function K1Q_X = getstartingvaluesformle(nSeeds,gamma,Y,W,cP,mats,dt,Sigma_cP)

N = size(cP,2);
LLK = zeros(nSeeds,1);
LAMQ = - sort( - log(unifrnd(0.9,1,N,nSeeds)));
parfor i = 1:nSeeds
    K1Q_X = diag(LAMQ(:,i));
    LLK(i) = jszllk(gamma,Y,W,cP,mats,dt,[],K1Q_X,Sigma_cP,[]);
end
K1Q_X = diag(LAMQ(:,LLK == min(LLK)));