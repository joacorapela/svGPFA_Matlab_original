function prs = Estep_Update_single_svGPFA(m,ntr);
% extract E step parameters to be optimized
mu0 = [];
sig0= [];
sigdiag0 = [];

for k = 1:m.dx
    mu0 = [mu0; m.q_mu{k}(:,:,ntr)];
    sig0 = [sig0;m.q_sqrt{k}(:,:,ntr)];
    sigdiag0 = [sigdiag0;m.q_diag{k}(:,:,ntr)];
end

% make objective function handle
prs0 = [mu0;sig0;sigdiag0];

% extract current hyperparameters
current_hprs = cellfun(@(struct)struct.hprs, m.kerns,'uni',0)'; % extract kernel hyperparams

% extract Z's for current trial
single_Z = cellfun(@(x) x(:,:,ntr),m.Z,'uni',0);

% build kernel matrices
Kmats = BuildKernelMatrices(m,m.tt,single_Z,current_hprs,0);

% make objective function
fun = @(prs) Estep_Objective_svGPFA(m,prs,Kmats,ntr);

% check gradients numerically
if m.opts.verbose == 2 % extra level of verbosity
    fprintf('E step Grad Check:\n')
    DerivCheck(fun,prs0);
end

% run optimizer
optimopts = optimset('Gradobj','on','display', 'none');
optimopts.MaxIter = m.opts.maxiter.Estep;

prs = minFunc(fun,prs0,optimopts);
