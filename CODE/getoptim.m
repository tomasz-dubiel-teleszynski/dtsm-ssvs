function pars = getoptim(pars,rest,satruefalse)

options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',1000,'MaxFunEvals',2000);

warning('off') % switch off

% optimization
for i = 1:3
    [pars,Llk0,~,~,~,hesspars] = fminunc(@objmle,pars,options,rest);
end

warning('on')  % switch on

% display likelihood
Llk0

if satruefalse  
    % ...    
end

end