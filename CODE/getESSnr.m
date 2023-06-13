function [ESSnr, Wnr] = getESSnr(fi,fi1,logiw)

logiw = logiw*(fi - fi1);
iw = exp(logiw);
Wnr = iw/sum(iw);
ESSnr = 1/sum(Wnr.^2);
% NaN and Inf checks
if isnan(ESSnr) || isinf(ESSnr)
    ESSnr = 1;
end

end