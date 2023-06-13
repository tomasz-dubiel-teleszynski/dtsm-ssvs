function A = makepd(A)
% This function makes matrix A semi-positive definite if it is not and also
% symmetric just in case

% check if semi-positive definite
[V,D] = eig(A);
d = diag(D);
if min(d) < 0
    % make semi-positive definite
    d(d < 0) = 1e-7;
    D = diag(d);
    A = V*D*V';
end

% make symmetric
A = (A + A')/2;

end