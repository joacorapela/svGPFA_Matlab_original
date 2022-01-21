function m = variationalEM(m);
% m = variationalEM(m);
%
% this function fits a sparse variational GPFA model to multivariate 
% observations using variational EM
% 
% input:
% ======
% m         -- model structure created using InitialiseModel_svGPFA (or
%              (other model initialiser)
%
% output:
% ======
% m         -- model structure with optimised parameters
%
% See also: InitialiseModel_svGPFA, InitialiseModel_grouped_svGPFA, 
%           InitialiseModel_warped_grouped_svGPFA
%
%
% Duncker, 2018
%
%
t_start = tic; % record starting time
abstol = 1e-05; % convergence tolerance
m.FreeEnergy = []; 

% print output
if m.opts.verbose
    fprintf('%3s\t%10s\t%10s\t%10s\n','iter','objective','increase','iterTime')
end

% ========= run variational inference =========
t_start_iter = tic;
for i = 1:m.opts.maxiter.EM
   
    % ========= resample minibatch ============    
    
    m = m.EMfunctions.sampleMinibatch(m);
    
    % ========= E-step: update variational parameters =========
    
    m = Estep(m);
    
    % ========= compute new value of free energy ==================
    
    m.FreeEnergy(i,1) = m.EMfunctions.VariationalFreeEnergy(m);
    m.iterationTime(i,1) = toc(t_start_iter);
    
    if i > 1
        FEdiff = m.FreeEnergy(i,1) - m.FreeEnergy(i-1,1);
    else
        FEdiff = NaN;
    end
    
    if m.opts.verbose
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\n',i,m.FreeEnergy(i),FEdiff,m.iterationTime(i));
    end
    
    % ========= check convergence in free energy =========
    if i > 2 && abs(m.FreeEnergy(i,1) - m.FreeEnergy(i-1,1)) < abstol
        break;
    end
    t_start_iter = tic;
    
    % ========= M-step: optimise wrt model parameters =========
    
    if i > 1 % skip first M-step to avoid early convergence to bad optima
        m = Mstep(m); 
    end

    % ========= hyper-M step: optimise wrt hyperparameters =========
    
    m = hyperMstep(m);
    
    % ========= inducing point hyper-M step: optimise wrt inducing point locations =========
    
    m = inducingPointMstep(m);
end

% save and report elapsed time
m.RunTime = toc(t_start);
if m.opts.verbose
    fprintf('Elapsed time is %1.5f seconds\n',m.RunTime);
end
