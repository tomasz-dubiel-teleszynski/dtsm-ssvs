function acr = acrate(theta)
% This function calculates acceptance rate for a matrix of posterior draws
% where each column relates to different parameters

% transpose if necessary
[rs,cs] = size(theta);
if rs < cs 
    theta = theta'; 
end

% get number of parameters
K = size(theta,2);

% calculate acceptance rate(s)
acr = zeros(1,K);
for i = 1:K
    k1 = (diff(theta(:,i)) ~= 0);
    acr(i) = mean(k1);
end

end