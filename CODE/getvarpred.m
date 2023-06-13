function cPpred = getvarpred(randnmat,h,W,cP,mats,dt,kinfQ,K1Q_X,lam0,lam1,Sigma_cP)

[K0Q_cP,K1Q_cP] = jszloadings(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);

N = size(cP,2);
cPpred = zeros(h,N);

[K0P_cP,I_K1P_cP] = getPdyn(K0Q_cP,K1Q_cP,lam0,lam1);

for i = 1:h
    if i > 1
        cPpred(i,:) = K0P_cP' + cPpred(i-1,:)*I_K1P_cP' + (Sigma_cP*randnmat(:,i))'; 
    else
        cPpred(i,:) = K0P_cP' + cP(end,:)*I_K1P_cP' + (Sigma_cP*randnmat(:,i))';  
    end
end

end