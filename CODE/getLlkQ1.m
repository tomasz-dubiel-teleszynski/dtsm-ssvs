function LlkQ = getLlkQ1(Y,W,cP,AcP,BcP,sige)
% This function calculates Q log for 1 observation

% get number of factors
N = size(cP,2);

% get number of observations and number of yields
[T,J] = size(Y);

% get Q log likelihood
Y_hat = repmat(AcP',T,1) + cP(end,:)*BcP';
%errors = (Y - Y_hat)*null(W); % null space of W
errors = Y - Y_hat;
llkQ = - (J - N)*0.5*log(2*pi) -(J - N)*log(sige) - 0.5*diag(errors*errors')/(sige*sige);
LlkQ = sum(llkQ);

end