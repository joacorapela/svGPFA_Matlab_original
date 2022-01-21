function prs = Estep_Update_all_PointProcess_svGPFA(m)

qmu0 = cell2mat(m.q_mu(:));
qsqrt0 = cell2mat(m.q_sqrt(:));
qdiag0 = cell2mat(m.q_diag(:));

prs0 = [qmu0;qsqrt0;qdiag0];
prs0 = prs0(:);

% extract current hyperparameters
current_hprs = cellfun(@(struct)struct.hprs, m.kerns,'uni',0)'; % extract kernel hyperparams

% get current Kernel matrices for quadrature 
Kmats_Quad = BuildKernelMatrices(m,m.ttQuad,m.Z,current_hprs,0);
% get current Kernel matrices evaluated at observed spike data
Kmats_Spikes = BuildKernelMatrices_fromSpikes(m,m.Z,current_hprs,0,1:m.ntr);

% make objective function
fun = @(prs) Estep_Objective_PointProcess_svGPFA(m,prs,Kmats_Quad,Kmats_Spikes);

% check gradients numerically
if m.opts.verbose == 2 % extra level of verbosity
    fprintf('E step Grad Check:\n')
    DerivCheck(fun,prs0);
end
% run optimizer
optimopts = optimset('Gradobj','on','display', 'none');
optimopts.MaxIter = m.opts.maxiter.Estep;

prs = minFunc(fun,prs0,optimopts);
