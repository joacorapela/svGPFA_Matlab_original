function prs = inducingPointMstep_all(m);

prs0 = vec(cell2mat(m.Z));

% make objective function
fun = @(prs) m.EMfunctions.InducingPointMstep_Objective(m,prs);
if m.opts.verbose == 2 % extra level of verbosity 
    fprintf('inducing point grad check:\n')
    DerivCheck(fun,prs0);
end
% run optimizer
optimopts = optimset('Gradobj','on','display', 'none');
optimopts.MaxIter = m.opts.maxiter.inducingPointMstep;

prs = minFunc(fun,prs0,optimopts);
