function [gamma, prob] = gammagsSSVS(lam0,lam1,...
                                     gamma,...
                                     EX_lam0lam1,...
                                     VAR_lam0lam1)

I = numel(gamma);
prob = zeros(1,I);
lam0lam1 = [lam0;lam1(:)];
for i = randperm(I)
    a = (1/VAR_lam0lam1.tau1(i))*exp(-0.5*(lam0lam1(i)^2)/(VAR_lam0lam1.tau1(i)^2))*EX_lam0lam1(i);
    b = (1/VAR_lam0lam1.tau0(i))*exp(-0.5*(lam0lam1(i)^2)/(VAR_lam0lam1.tau0(i)^2))*(1 - EX_lam0lam1(i));
    prob(i) = a/(a + b);
    gamma(i) = (rand() < prob(i));
end
