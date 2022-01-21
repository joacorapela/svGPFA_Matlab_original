function [obj,grad] = Estep_Objective_svGPFA(m,prs,Kmats,ntr);

if nargin < 4
    trEval = 1:m.ntr;
elseif length(ntr) == 1
    trEval = ntr;
else
    trEval = ntr;
end

[q_mu, q_sqrt, q_diag, q_sigma, idx, idx_sig,idx_sigdiag] = extract_variationalPrs_svGPFA(m,prs,trEval);

[mu_h,var_h] = m.EMfunctions.predict_MultiOutputGP(Kmats,q_mu,q_sigma,m.prs.C,m.prs.b);

% get expected log-likelihood and gradient
Elik = m.EMfunctions.likelihood(m,mu_h,var_h,trEval);
gradElik = m.EMfunctions.gradLik_variationalPrs(m,Kmats,q_sqrt, q_diag,idx,idx_sig,idx_sigdiag,trEval,mu_h,var_h);

% get KL divergence and gradient
KLd = build_KL_divergence(m,Kmats,q_mu,q_sigma);
gradKLd = grad_variationalPrs_KL_divergence(m,Kmats,q_mu,q_sqrt,q_diag,idx,idx_sig,idx_sigdiag,trEval);

obj = -Elik + KLd; % negative free energy
grad = -gradElik + gradKLd; % gradients
