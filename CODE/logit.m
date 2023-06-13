function res = logit(x,Min,Max)

if nargin == 1
    Min = 0*ones(size(x));
    Max = 1*ones(size(x));
end

res = min(max(-10^100,log((x-Min)./(Max-x))),10^100);
