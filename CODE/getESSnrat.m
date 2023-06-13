function [ESSnr, Wnr] = getESSnrat(fi,fi1,logiw,logwnr)

logwnr = logwnr + logiw*(fi - fi1);
wnr = exp(logwnr);
Wnr = wnr/sum(wnr);
ESSnr = 1/sum(Wnr.^2);
% NaN and Inf checks
if isnan(ESSnr) || isinf(ESSnr)
    ESSnr = 1;
end

end