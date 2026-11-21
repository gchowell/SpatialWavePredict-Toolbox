function [Phatss,diagnostics] = spatialWaveBootstrapRefits(context,curves,P,LB,UB,isLog,startSets,settings,useParallel)
%SPATIALWAVEBOOTSTRAPREFITS Refit every bootstrap dataset of one ranked model.
%   [Phatss,diagnostics] = spatialWaveBootstrapRefits(context,curves,P,LB,UB, ...
%       isLog,startSets,settings,useParallel)
%
%   context    fit context of the ranked model (spatialWaveObjective.m)
%   curves     the synthetic datasets, one per column
%   P          the ranked fit, [r p a K q alpha d]
%   LB, UB, isLog   search space (spatialWaveSearchSpace.m)
%   startSets  cell array, one matrix of start points per dataset (rows in
%              original units), drawn beforehand on the MATLAB client so that
%              the result does not depend on how the work is shared out
%   settings   spatialWaveFitSettings.m; each refit makes one MultiStart pass
%              from its own start points, serially
%   useParallel  true: the datasets are refitted on the parallel pool
%
%   Phatss       one refit per row
%   diagnostics  struct of column vectors, one entry per dataset: success,
%                usedFallback, exitflag, nActive, startAgreement,
%                gradientInfNorm, polishRun, polishImproved, nAtBound, fval
%
%   RETENTION POLICY (as before). A refit with no converged local minimum is
%   counted as unresolved and its best valid end point is kept. A dataset
%   with no valid parameter vector at all stops the run.
%
%   The loop over the datasets is the only parfor of the toolbox. It is in
%   this file so that it runs on the client, as an ordinary loop, whenever
%   useParallel is false.
%
%   Revision SW-2.
%
%   See also fittingModifiedLogisticFunctionPatchMultiple, fitSpatialWaveModel

M = size(curves,2);
refitSettings = settings;
refitSettings.useParallel = false;   % no pool inside a pool
refitSettings.fitMaxAttempts = 1;    % start points are given; none are drawn on a worker

Phatss = zeros(M,7);
success = false(M,1);
usedFallback = false(M,1);
exitflag = NaN(M,1);
nActive = NaN(M,1);
startAgreement = NaN(M,1);
gradientInfNorm = NaN(M,1);
polishRun = false(M,1);
polishImproved = false(M,1);
nAtBound = NaN(M,1);
fvals = NaN(M,1);
failed = false(M,1);
messages = cell(M,1);

nWorkers = 0;
if useParallel
    try
        pool = gcp();
        nWorkers = pool.NumWorkers;
    catch poolError
        warning('SpatialWave:ParallelUnavailable', ...
            'Parallel pool unavailable (%s); the bootstrap refits run serially.',poolError.message);
        nWorkers = 0;
    end
end

parfor (b = 1:M, nWorkers)
    refitContext = context;
    refitContext.ydata = curves(:,b);
    [Pb,fb,info] = fitSpatialWaveModel(refitContext,P,LB,UB,isLog,startSets{b},0,refitSettings);
    nActiveB = info.nActive;
    fallbackB = false;
    failedB = false;
    messageB = '';
    if ~info.success
        if ~isempty(info.fallback)
            % unresolved refit: keep its best valid end point
            Pb = info.fallback.P;
            fb = info.fallback.fval;
            nActiveB = info.fallback.nActive;
            fallbackB = true;
        else
            failedB = true;
            messageB = info.message;
            Pb = NaN(1,7);
        end
    end
    Phatss(b,:) = Pb;
    success(b) = info.success;
    usedFallback(b) = fallbackB;
    exitflag(b) = info.exitflag;
    nActive(b) = nActiveB;
    startAgreement(b) = info.startAgreement;
    gradientInfNorm(b) = info.gradientInfNorm;
    polishRun(b) = info.polish.run;
    polishImproved(b) = info.polish.improved;
    nAtBound(b) = numel(info.atLowerBound)+numel(info.atUpperBound);
    fvals(b) = fb;
    failed(b) = failedB;
    messages{b} = messageB;
end

if any(failed)
    first = find(failed,1);
    error('SpatialWave:BootstrapFitFailed', ...
        'Bootstrap dataset %d has no valid parameter vector: %s',first,messages{first});
end

diagnostics = struct('success',success,'usedFallback',usedFallback,'exitflag',exitflag, ...
    'nActive',nActive,'startAgreement',startAgreement,'gradientInfNorm',gradientInfNorm, ...
    'polishRun',polishRun,'polishImproved',polishImproved,'nAtBound',nAtBound,'fval',fvals);
end
