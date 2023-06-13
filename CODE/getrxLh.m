function rxL_ = getrxLh(h,AcPL,BcPL,YL_T_h,cP_T_h,mats,matsL)

% get number of observations
T_h = size(YL_T_h,1);

% construct matrix of maturities 
n = mats(find(mats==24):end);

n_h = n - h;

n = repmat(n,numel(h),1);
n_h = n_h(:);

h = repmat(h,numel(mats(find(mats==24):end)),1);
h  = h(:);

% number of excess returns
K = numel(n);

% calculate excess returns and create their labels
rx = zeros(T_h,K);
rxa = rx;
rxl = cell(K,1);

% containers for holding-periods
rxh = zeros(K,1);
rxn_h = rxh;
rxn = rxh;

for k = 1:K
    
    % auxiliary expressions
    tmp1    = h(k);
    tmp1idx = ismember(matsL,tmp1);
    tmp2    = n_h(k);
    tmp2idx = ismember(matsL,tmp2);
    tmp3    = n(k);
    tmp3idx = ismember(matsL,tmp3);
    
    shift = max(h) - tmp1;
    Yk  = YL_T_h(1:end-shift,:);
    cPk = cP_T_h(1:end-shift,:);
    
    % basic formula
    rx(:,k) = [NaN(tmp1,1); - tmp2*Yk(tmp1+1:end,tmp2idx) + tmp3*Yk(1:end-tmp1,tmp3idx) ...
        - tmp1*Yk(1:end-tmp1,tmp1idx);NaN(shift,1)];

    % alternative formula
    rxa(:,k) = [NaN(tmp1,1);(AcPL(tmp2idx)*(-tmp2) - AcPL(tmp3idx)*(-tmp3) + AcPL(tmp1idx)*(-tmp1) + ...
        BcPL(tmp2idx,:)*(-tmp2)*cPk(tmp1+1:end,:)' - ...
        (BcPL(tmp3idx,:)*(-tmp3) - BcPL(tmp1idx,:)*(-tmp1))*cPk(1:end-tmp1,:)')';NaN(shift,1)];
    
    % store holding-periods
    rxh(k)   = tmp1;
    rxn_h(k) = tmp2;
    rxn(k)   = tmp3;
    
    % common labels
    rxl{k} = ['h-',num2str(tmp1),'-n-h-',num2str(tmp2),'-n-',num2str(tmp3)];

end

% means
rxbar = zeros(1,K);
rxabar = zeros(1,K);
for k = 1:K
    tmp = h(k);
    shift = max(h) - tmp;
    rxbar(k)  = mean(rx(h(k)+1:end-h(k)-shift,k),'omitnan');
    rxabar(k) = mean(rxa(h(k)+1:end-h(k)-shift,k),'omitnan');   
end

% run predictive regressions
R2 = zeros(K,1);
for k = 1:K
    tmp = h(k);
    shift = max(h) - tmp;
    cPk = cP_T_h(1:end-shift,:);
    [~,~,~,~,stats] = regress(rx(tmp+1:end-shift,k),[ones(T_h-tmp-shift,1),cPk(1:end-tmp,:)],0.05);
    R2(k) = stats(1);
end

% allocate forecasts to the end
for k = 1:K
    tmp = h(k);
    shift = max(h) - tmp;
    rx(end,k)  = rx(end-shift,k);
    rxa(end,k) = rxa(end-shift,k); 
end

% load rx_ structure
rxL_          = struct();
rxL_.rx       = rx;
rxL_.rxbar    = rxbar;
rxL_.rxa      = rxa;
rxL_.rxabar   = rxabar;
rxL_.rxl      = rxl;
rxL_.rxh      = rxh;
rxL_.rxn_h    = rxn_h;
rxL_.rxn      = rxn;
rxL_.R2       = R2;