function logProp = mvtlogProp(X,mu,cov_mat,df)
% This function calculates log proposal density as log MVT density

rescale = sqrt(diag(cov_mat));
corr_mat = cov_mat ./ (rescale * rescale');

X = (X - mu) ./ rescale;

logProp = log(mvtpdf(X,corr_mat,df));

end