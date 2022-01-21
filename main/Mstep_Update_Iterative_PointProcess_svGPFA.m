function m = Mstep_Update_Iterative_PointProcess_svGPFA(m)

% get current parameters
prs0 = [m.prs.C(:);m.prs.b];

% build kernel matrices
current_hprs = cellfun(@(struct)struct.hprs, m.kerns,'uni',0)'; % extract kernel hyperparams
Kmats_Quad = BuildKernelMatrices(m,m.ttQuad,m.Z,current_hprs,0); % get current Kernel matrices  
Kmats_Spikes = BuildKernelMatrices_fromSpikes(m,m.Z,current_hprs,0); % get current Kernel matrices evaluated at observed spike data

% predict posterior means and variances
[mu_k_Quad, var_k_Quad] = predict_posteriorGP(Kmats_Quad,m.q_mu,m.q_sigma);
[mu_k_Spikes, var_k_Spikes] = predict_posteriorGP_fromSpikes(Kmats_Spikes,Kmats_Quad.Kzzi,m.q_mu,m.q_sigma,1:m.ntr);

% make objective function
fun = @(prs) Mstep_Objective_PointProcess(m,prs,mu_k_Quad,var_k_Quad,mu_k_Spikes,var_k_Spikes);

% check gradients numerically
if m.opts.verbose == 2 % extra level of verbosity
    fprintf('M step Grad Check:\n')
    DerivCheck(fun,prs0);
end
% minimize
optimopts = optimset('Gradobj','on','display', 'none');
optimopts.MaxIter = m.opts.maxiter.Mstep;
prs = minFunc(fun,prs0,optimopts);

% update model parameters in structure
m.prs.C = reshape(prs(1:m.dy*m.dx),[m.dy, m.dx]);
m.prs.b = prs(m.dy*m.dx + 1 : end);