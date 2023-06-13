function gamma = gammarest(model,N)

% restrictions on lam0 and lam1
gamma = zeros(1,N*(N+1));
if N == 3
    switch model
        case 'LAM0_1_LAM1_11_LAM1_12_LAM1_23'
            gamma(1)  = 1;  % lam0(1)
            gamma(4)  = 1;  % lam1(1,1)
            gamma(7)  = 1;  % lam1(1,2)
            gamma(11) = 1;  % lam1(2,3)
        case 'LAM0_1_LAM1_11_LAM1_12'
            gamma(1)  = 1;  % lam0(1)
            gamma(4)  = 1;  % lam1(1,1)
            gamma(7)  = 1;  % lam1(1,2)
        case 'LAM0_1_LAM1_12'
            gamma(1)  = 1;  % lam0(1)
            gamma(7)  = 1;  % lam1(1,2)
        case 'LAM1_11_LAM1_12'
            gamma(4)  = 1;  % lam1(1,1)
            gamma(7)  = 1;  % lam1(1,2)
        case 'LAM1_11'
            gamma(4)  = 1;  % lam1(1,1)
        case 'LAM1_12'
            gamma(7)  = 1;  % lam1(1,2)
        case 'MAXIMALLY_FLEXIBLE'
            gamma = ones(1,N*(N+1));
        case 'LAM1_13_LAM1_23_ZEROS'
            gamma = ones(1,N*(N+1));
            gamma(10) = 0;  % lam1(1,3)
            gamma(11) = 0;  % lam1(2,3)
        case 'CONSTANT_TERM_PREMIA'
            gamma(1)  = 1;  % lam0(1)
            gamma(2)  = 1;  % lam0(2)
            gamma(3)  = 1;  % lam0(3)
        case 'ZERO_TERM_PREMIA'
            
    end
    % % non-zero restrictions
    % % lam0
    % gamma(1)  = 0;  % lam0(1)
    % gamma(2)  = 0;  % lam0(2)
    % gamma(3)  = 0;  % lam0(3)
    % % lam1 (columnwise)
    % gamma(4)  = 0;  % lam1(1,1)
    % gamma(5)  = 0;  % lam1(2,1)
    % gamma(6)  = 0;  % lam1(3,1)
    % gamma(7)  = 0;  % lam1(1,2)
    % gamma(8)  = 0;  % lam1(2,2)
    % gamma(9)  = 0;  % lam1(3,2)
    % gamma(10) = 0;  % lam1(1,3)
    % gamma(11) = 0;  % lam1(2,3)
    % gamma(12) = 0;  % lam1(3,3)
elseif N == 2
    switch model
        case 'LAM0_1_LAM1_11_LAM1_12'
            gamma(1)  = 1;  % lam0(1)
            gamma(3)  = 1;  % lam1(1,1)
            gamma(5)  = 1;  % lam1(1,2)
        case 'LAM0_1_LAM1_12'
            gamma(1)  = 1;  % lam0(1)
            gamma(5)  = 1;  % lam1(1,2)
        case 'LAM1_11_LAM1_12'
            gamma(3)  = 1;  % lam1(1,1)
            gamma(5)  = 1;  % lam1(1,2)
        case 'LAM1_11'
            gamma(3)  = 1;  % lam1(1,1)
        case 'LAM1_12'
            gamma(5)  = 1;  % lam1(1,2)
        case 'MAXIMALLY_FLEXIBLE'
            gamma = ones(1,N*(N+1));
    end
    % % non-zero restrictions
    % % lam0
    % gamma(1)  = 0;  % lam0(1)
    % gamma(2)  = 0;  % lam0(2)
    % % lam1 (columnwise)
    % gamma(3)  = 0;  % lam1(1,1)
    % gamma(4)  = 0;  % lam1(2,1)
    % gamma(5)  = 0;  % lam1(1,2)
    % gamma(6)  = 0;  % lam1(2,2)
end

end