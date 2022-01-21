function [obj,grad] = Estep_Objective_PointProcess_svGPFA(m,prs,Kmats_Quad,Kmats_Spikes,ntr);

if nargin < 5
    trEval = 1:m.ntr;
elseif length(ntr) == 1
    trEval = ntr;
else
    trEval = ntr;
end

[q_mu, q_sqrt, q_diag, q_sigma, idx, idx_sig,idx_sigdiag] = extract_variationalPrs_svGPFA(m,prs,trEval);

% get multi output prediction for quad and from spikes
[mu_h_Quad,var_h_Quad] = m.EMfunctions.predict_MultiOutputGP(Kmats_Quad,q_mu,q_sigma,m.prs.C,m.prs.b);
[mu_h_Spikes,var_h_Spikes] = m.EMfunctions.predict_MultiOutputGP_fromSpikes(Kmats_Spikes,Kmats_Quad.Kzzi,q_mu,q_sigma,m.prs.C,m.prs.b,trEval,m.index);

% get expected log-likelihood and gradient
Elik = m.EMfunctions.likelihood(m,mu_h_Quad,var_h_Quad,mu_h_Spikes,var_h_Spikes,trEval);
gradElik = m.EMfunctions.gradLik_variationalPrs(m,Kmats_Quad,Kmats_Spikes,q_sqrt, q_diag,idx,idx_sig,idx_sigdiag,trEval,...
    mu_h_Quad,var_h_Quad);

% get KL divergence and gradient
KLd = build_KL_divergence(m,Kmats_Quad,q_mu,q_sigma);
gradKLd = grad_variationalPrs_KL_divergence(m,Kmats_Quad,q_mu,q_sqrt,q_diag,idx,idx_sig,idx_sigdiag,trEval);

% assemble objective function and gradients 
obj = -Elik + KLd; % negative free energy
grad = -gradElik + gradKLd; % gradients

