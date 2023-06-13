function [W,Y,cP,dates] = realWYcP(mats,N,varargin)
% This function calculated real factor loading matrix, real yields and
% observed factors from real yields

% dummies
Y = [];
dates = [];

% load real yields
load('yields_data');

% annualize
n_per = 12;
Y = Y(:,mats)/n_per;

% constrain time frame as in Bauer (2016)
% or provide own time span
if isempty(varargin)
    start_sample = 19900101;
    end_sample = 20080000;
else
    start_sample = varargin{1};
    end_sample = varargin{2};
end
sel_sample =  (dates >= start_sample) & (dates <= end_sample);
Y = Y(sel_sample,:);
dates = dates(sel_sample,:);

% do PCA with normalizations
[W,cP] = getWcP(Y,N);

end