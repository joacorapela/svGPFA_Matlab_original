% mu_k and var_k are the mean and covariance of the sparce posterior GP: 
% Eq. (5) from Dunker and Sahani, 2008 (completes de inner terms computed by predict_MultiOutputGP). Also see, Eq.6 from Titsias, 2009.

function [mu_h,var_h,mu_k,var_k] = predict_MultiOutputGP(Kmats,q_mu,q_sigma,C,b)

[mu_k,var_k] = predict_posteriorGP(Kmats,q_mu,q_sigma);

mu_h = bsxfun(@plus,mtimesx(mu_k, C'), b');
var_h = mtimesx(var_k,(C.^2)');
