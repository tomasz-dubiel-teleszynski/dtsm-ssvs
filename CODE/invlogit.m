function res = invlogit(x,Min,Max)

if nargin == 1
    Min = 0*ones(size(x));
    Max = 1*ones(size(x));
end

res = (Max.*exp(x)+Min)./(1+exp(x));
inds = find(or(isnan(res),isinf(res)));
res(inds) = (Max*ones(size(inds))+Min*ones(size(inds)).*exp(-x(inds)))./(1+exp(-x(inds)));

