
%% script to simulate Poisson Process data and run sv-ppGPFA 
% clear all; close all;
olddir = cd('..');
set_paths;
cd(olddir);

emMaxIter=-1; pEstNumber=86915600;
pythonDataFilenamePattern = '../data/%08d_estimationDataForMatlab.mat';

pythonDataFilename = sprintf(pythonDataFilenamePattern, pEstNumber);
pythonData = load(pythonDataFilename);

if(emMaxIter>=0)
    pythonData.emMaxIter = emMaxIter;
end

% control variables
dx = pythonData.nLatents;
dy = pythonData.nNeurons;
ntr = pythonData.nTrials;

% trials lengths
trLen = pythonData.trialsLengths;

% spikes
Y = cell(ntr,1);
for r=0:ntr-1
    Y{r+1,1} = cell(dy,1);
    for n=0:dy-1
        spikes = pythonData.(sprintf('spikesTimes_%d_%d', r, n));
        Y{r+1,1}{n+1,1} = reshape(spikes, length(spikes), 1);
    end
end

indPointsLocsGramMatrixEpsilon=pythonData.indPointsLocsKMSRegEpsilon;

% embedding params
noisyPRS.C = pythonData.C;
noisyPRS.b = pythonData.d;

% inducing points
Z = cell(dx,1);
for ii = 0:dx-1
    Zii = pythonData.(sprintf('indPointsLocs_%d', ii));
    Z{ii+1} = zeros(size(Zii,2),1,ntr);
    for jj = 0:ntr-1
        Z{ii+1}(:,1,jj+1) = Zii(jj+1,:);
    end
end

% kernels
kerns = {};
for ii = 0:dx-1
    kernelType = pythonData.(sprintf('kernelType_%d', ii));
    kernelsParams0 = pythonData.(sprintf('kernelsParams_%d', ii))';
    switch kernelType
        case 'PeriodicKernel'
            kerns{ii+1} = buildKernel('Periodic', kernelsParams0);
        case 'ExponentialQuadraticKernel'
            kerns{ii+1} = buildKernel('RBF', kernelsParams0);
        otherwise
            error(sprintf('kernelType %s not recognized', kernelType))
    end
end

%% initialise model structure
options.parallel = 0;
options.verbose = 1;

options.maxiter.EM = pythonData.emMaxIter;
options.maxiter.Estep = pythonData.eStepMaxIter;
options.maxiter.Mstep = pythonData.mStepEmbeddingMaxIter;
options.maxiter.hyperMstep = pythonData.mStepKernelsMaxIter;
options.maxiter.inducingPointMstep = pythonData.mStepIndPointsMaxIter;

options.nbatch = ntr;
options.nquad = length(pythonData.('legQuadPoints'));

m = InitialiseModel_svGPFA('PointProcess',@exponential,Y,trLen,kerns,Z,noisyPRS,options);
m.epsilon = indPointsLocsGramMatrixEpsilon; % value of diagonal added to kernel inversion for stability

% q_mu, q_sqrt, q_diag
q_mu = cell(1,dx);
q_sqrt = cell(1,dx);
q_diag = cell(1,dx);
for ii=0:dx-1
    nIndPointsk = size(Z{ii+1}, 1);
    q_mu{ii+1} = zeros(nIndPointsk, 1, ntr);
    q_mu{ii+1}(:,1,:) = pythonData.(sprintf('qMu_%d', ii))';
    q_sqrt{ii+1} = zeros(nIndPointsk, 1, ntr);
    q_sqrt{ii+1}(:,1,:) = pythonData.(sprintf('qSVec_%d', ii))';
    q_diag{ii+1} = zeros(nIndPointsk, 1, ntr);
    q_diag{ii+1}(:,1,:) = pythonData.(sprintf('qSDiag_%d', ii))';
end
q_sigma = get_full_from_lowplusdiag(m,q_sqrt,q_diag);
m.q_mu = q_mu;
m.q_sqrt = q_sqrt;
m.q_diag = q_diag;
m.q_sigma = q_sigma;

% ttQuad, wwQuad
ttQuad = zeros(length(pythonData.('legQuadPoints')), 1, ntr);
ttQuad(:,1,:) = pythonData.('legQuadPoints')';
wwQuad(:,1,:) = pythonData.('legQuadWeights')';
m.ttQuad = ttQuad;
m.wwQuad = wwQuad;

%% set extra options and fit model
% m.opts.maxiter.EM = 50; % maximum number of iterations to run
m.opts.fixed.Z = 0; % set to 1 to hold certain parameters values fixed
m.opts.fixed.hprs = 0;
m.opts.nbatch = ntr; % number of trials to use for hyperparameter update

m = variationalEM(m);

%% predict latents and MultiOutput GP
latentsTimes = pythonData.(sprintf('latentsTrialsTimes_0'));
latentsTimes = reshape(latentsTimes, length(latentsTimes), 1);
pred = predictNew_svGPFA(m, latentsTimes);

lowerBound = m.FreeEnergy;
% elapsedTime = m.elapsedTime;
meanEstimatedLatents = pred.latents.mean;
varEstimatedLatents = pred.latents.variance;

estimationResFilename = sprintf('results/%08d-pointProcessEstimationRes.mat', estResNumber);
save(estimationResFilename, 'lowerBound', 'latentsTimes', 'meanEstimatedLatents', 'varEstimatedLatents', 'm');
