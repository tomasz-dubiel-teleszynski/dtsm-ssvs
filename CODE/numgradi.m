function df = numgradi(func, X, h,varargin)
% This function calculates numerical gradient df for function func at X
% with tolerance h and using inputs to func in varargin

% get number of parameters
L = length(X);

% set container
df = zeros(L,1);

% for each dimension of objective function
for i=1:length(X)
    % vary variable i by a small amount (left and right)
    x1 = X;
    x2 = X;
    x1(i) = X(i) - h(i);
    x2(i) = X(i) + h(i);
    
    % evaluate the objective function at the left and right points
    y1 = feval(func,x1,varargin{:});
    y2 = feval(func,x2,varargin{:});
    
    % calculate the slope (rise/run) for dimension i
    df(i) = (y2 - y1)/(2*h(i));
end

end