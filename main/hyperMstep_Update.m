function m = hyperMstep_Update(m)
% function to update hyperparameters of svGPFA model

% extract hyperparameters for each GP kernel
prs0 = m.EMfunctions.extract_hyperParams(m);

% make objective function
fun = @(prs,nIter) m.EMfunctions.HyperMstep_Objective(m,prs,nIter);

% check gradients numerically
if m.opts.verbose == 2 % extra level of verbosity
    fprintf('hyper M step grad check:\n')
    DerivCheck(@(x)fun(x,1),prs0);
end

optimopts = optimset('Gradobj','on','display', 'none');

optimopts.MaxIter = m.opts.maxiter.hyperMstep;

% ADAM option values
DEF_stepSize = 0.001;
DEF_beta1 = 0.9;
DEF_beta2 = 0.999;
DEF_epsilon = sqrt(eps);

prs = fmin_adam(fun,prs0,DEF_stepSize, DEF_beta1, DEF_beta2, DEF_epsilon, 1, optimopts);

% update kernel hyperparameters in model structure
m = m.EMfunctions.updateHyperParams(m,prs);