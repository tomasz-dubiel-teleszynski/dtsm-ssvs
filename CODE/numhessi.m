function H = numhessi(func, X, h,varargin)
% This function calculates numerical hessian H for function func at X with
% tolerance h and using inputs to func in varargin

% get number of parameters
L = length(X);

% set container
H = zeros(L);

% for each dimension of objective function
for i=1:L
    % derivative at first point (left)
    x1 = X;
    x1(i) = X(i) - h(i);
    df1 = numgradi(func, x1, h, varargin{:});
    
    % derivative at second point (right)
    x2 = X;
    x2(i) = X(i) + h(i);
    df2 = numgradi(func, x2, h, varargin{:});
    
    % differentiate between the two derivatives
    d2f = (df2 - df1)/(2*h(i));
    
    % assign as row i of Hessian
    H(i,:) = d2f';
end

end