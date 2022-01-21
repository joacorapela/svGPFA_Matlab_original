function pred = predictNew_svGPFA(m,testTimes);

current_hprs = cellfun(@(struct)struct.hprs, m.kerns,'uni',0)'; % extract kernel hyperparams
Kmats = BuildKernelMatrices(m,testTimes,m.Z,current_hprs,0); % get current Kernel matrices   

[mu_k, var_k] = predict_posteriorGP(Kmats,m.q_mu,m.q_sigma);

mu_h = bsxfun(@plus,mtimesx(mu_k,m.prs.C'), m.prs.b');
var_h = mtimesx(var_k,(m.prs.C.^2)');

pred.latents.mean = mu_k;
pred.latents.variance = var_k;
pred.multiOutputGP.mean = mu_h;
pred.multiOutputGP.variance = var_h;