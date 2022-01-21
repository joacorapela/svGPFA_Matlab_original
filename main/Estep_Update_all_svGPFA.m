function prs = Estep_Update_all_svGPFA(m)

qmu0 = cell2mat(m.q_mu(:));
qsqrt0 = cell2mat(m.q_sqrt(:));
qdiag0 = cell2mat(m.q_diag(:));

prs0 = [qmu0;qsqrt0;qdiag0];
prs0 = prs0(:);

% extract current hyperparameters
current_hprs = cellfun(@(struct)struct.hprs, m.kerns,'uni',0)'; % extract kernel hyperparams

% build kernel matrices
Kmats = BuildKernelMatrices(m,m.tt,m.Z,current_hprs,0);

% make objective function
fun = @(prs) Estep_Objective_svGPFA(m,prs,Kmats);

% check gradients numerically
if m.opts.verbose == 2 % extra level of verbosity
    fprintf('E step Grad Check:\n')
    DerivCheck(fun,prs0);
end
% run optimizer
optimopts = optimset('Gradobj','on','display', 'none');
optimopts.MaxIter = m.opts.maxiter.inducingPointMstep;

prs = minFunc(fun,prs0,optimopts);
