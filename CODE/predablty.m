function [rxf,rxo,rxobar,rxl,rxfall,rxoall] = predablty(t,ws1,thetast1,dims,YL,matsL,W,cP,mats,dt)

% number of factors
N = size(cP,2);

% horizon for VAR predictions
h = 1:12;

% choose relevant observations
cPt    = cP(1:t,:); 
YL_t_h = YL(1:t+h(end),:);
cP_t_h = cP(1:t+h(end),:);
cP_t_ho = cP_t_h;

Np = size(thetast1,1);

rxf = zeros(Np,(numel(mats(find(mats==12):end))-1)*numel(h));

rxo = rxf;
rxobar = rxf;
rxl = cell(Np,1);

rxoall = cell(Np,1);

% IN HERE THERE IS A RANDOM NUMBER GENERATION STEP
% SO THAT WITH FIXED SEED PREDICTIONS DON'T CHANGE
randnmat = cell(Np,1);
for ip = 1:Np
   randnmat{ip} = randn(N,h(end)); 
end

parfor ip = 1:Np    
    
    % transform back parameters from single particle
    [kinfQ,K1Q_X,lam0,lam1,Sigma_cP] = pars2svms(thetast1(ip,:)',dims);
    kinfQ = rescalekinfQ(kinfQ);
    K1Q_X = mlgitinvs(K1Q_X);
    Sigma_cP = expd(Sigma_cP);
    
    % get AcPL, BcPL
    [~,~,~,~,~,~,~,~,AcPL,BcPL] = jszloadingspost(W,mats,dt,kinfQ,K1Q_X,Sigma_cP);

    if ip == 1
        % get observed excess returns
        rxoLh_ = getrxLh(h,AcPL,BcPL,YL_t_h,cP_t_ho,mats,matsL);
        rxo(ip,:) = rxoLh_.rx(end,:);
        rxoall{ip} = rxoLh_.rx(2:end-12,:);
        rxobar(ip,:) = rxoLh_.rxbar;
        rxl{ip} = rxoLh_.rxl;
    end
    
    % get VAR predictions
    cPpred = getvarpred(randnmat{ip},h(end),W,cPt,mats,dt,kinfQ,K1Q_X,lam0,lam1,Sigma_cP);
    cP_t_h = [cPt;cPpred];
    
    % get implied excess returns
    rxfLh_ = getrxLh(h,AcPL,BcPL,YL_t_h,cP_t_h,mats,matsL);
    rxf(ip,:) = rxfLh_.rxa(end,:);
 
end

rxl = rxl{1};
rxo = rxo(1,:);
rxobar = rxobar(1,:);
rxfall = rxf;
rxoall = rxoall{1};
rxf = sum(rxf.*repmat(ws1,1,numel(rxo)))/sum(ws1);

end