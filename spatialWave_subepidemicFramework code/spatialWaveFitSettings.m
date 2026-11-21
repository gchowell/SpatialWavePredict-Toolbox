function S = spatialWaveFitSettings()
%SPATIALWAVEFITSETTINGS Optimizer and search-space settings of the fitting layer.
%   S = spatialWaveFitSettings() returns the defaults below. If a function
%   options_fit_settings_sw.m is on the path, it is called with no inputs and
%   must return a struct; each of its fields replaces the default of the same
%   name. Unknown field names are an error, so a misspelled setting cannot be
%   ignored silently. Copy options_fit_settings_sw_example.m to
%   options_fit_settings_sw.m to start. Existing options*.m files need no
%   change. (The name ends in _sw so that it cannot be mistaken for the
%   settings file of the SubEpiPredict toolbox when both are on the path.)
%
%   SEARCH SPACE (spatialWaveSearchSpace.m)
%     logScaleR, logScaleK, logScaleAlpha   true: r, K and the negative-
%         binomial alpha are searched on a log scale. The bounds do not
%         change; only the coordinates the optimizer and the start points
%         work in.
%     rLB, rUB, pLB, pUB, aLB, aUB, KLB, KUB, qLB, qUB   [] = the earlier
%         bound (getSpatialWaveBounds.m); a number replaces it. A lower bound
%         equal to its upper bound fixes the parameter.
%     rUpperBoundFactor   10: the upper bound of r is at least 10 times the
%         largest count of the calibration series. 0: the earlier bound alone.
%
%   FITTED SERIES (spatialWaveFitSeries.m)
%     countLikelihoodOnRawData   true: Poisson and negative-binomial fits
%         (methods 1, 3, 4, 5) use the raw counts whatever smoothfactor1 is;
%         smoothfactor1 then only shapes the grid of onset thresholds.
%         false: they are fitted to the moving average, as before.
%
%   OPTIMIZER (fitSpatialWaveModel.m)
%     optimOptions      struct of fmincon option name/value pairs applied on
%                       top of the defaults (sqp, StepTolerance 1e-10,
%                       FunctionTolerance 1e-10, OptimalityTolerance 1e-8,
%                       forward differences with step 1e-6 in unit coordinates).
%     fitMaxAttempts    3: MultiStart passes with fresh start points before a
%                       candidate is reported as failed.
%     useParallel       false. true: the local solves of each MultiStart call
%                       of the search, and the bootstrap refits, run on a
%                       parallel pool (Parallel Computing Toolbox).
%     polish            true: Nelder-Mead pass after MultiStart, run only when
%                       the stationarity check fails.
%     polishGradientTol 10: the polish runs when the largest gradient
%                       component, in log-likelihood units per unit
%                       coordinate, exceeds this value.
%     polishBudget      100: function evaluations of the polish per free
%                       parameter.
%     agreementTol      1e-2: two local runs agree when their objectives
%                       differ by less than this many log-likelihood units.
%
%   SEARCH (fittingModifiedLogisticFunctionPatchABC.m)
%     continuationStart true: the fit of each onset threshold also starts
%                       from the accepted solution of the previous threshold.
%     firstThresholdStartFactor 3: a fit with no solution of a previous
%                       threshold to start from gets this many times
%                       numstartpoints start points. 1: no extra start points.
%     thresholdLowPoints 6: onset thresholds added, evenly on a log scale,
%                       between the first observation and the first step of
%                       the regular grid. The regular grid is evenly spaced in
%                       cumulative cases, so its first step can be thousands
%                       of cases wide, and waves whose sub-epidemics start
%                       early were not among the candidates. 0: the regular
%                       grid alone, as before.
%
%   BOOTSTRAP (fittingModifiedLogisticFunctionPatchMultiple.m)
%     bootstrapStarts   1: Latin-hypercube start points of each bootstrap
%                       refit, besides the ranked fit and one jittered copy
%                       of it.
%
%   REPEATABLE RUNS (spatialWaveSeed.m)
%     rngSeed           1: the random number generator is reset with this
%                       seed at the start of the fit (once per area) and of
%                       the forecast script, so the same call on the same
%                       files gives the same numbers. [] leaves the generator
%                       alone, as before revision SW-2.
%
%   FORECAST SCORES (plotForecast_SW_subepidemicFramework.m)
%     baselineWindow    []: the flat baseline (spatialWaveBaselineForecast.m)
%                       repeats the mean of the last 7 observations of daily
%                       data and the last observation of weekly or yearly
%                       data. A positive integer replaces that number. It is
%                       not a fit setting: changing it does not make a saved
%                       fit out of date. (Revision SW-3.)
%
%   IDENTIFIABILITY REPORT (spatialWaveIdentifiability.m; revision SW-6). These
%   settings change the report only, not a fit: saved fits stay valid.
%     PIIscale          []: the headline measure of each parameter is chosen
%                       by its lower search bound (positive: the PII on the
%                       log scale; otherwise the share of the search range).
%                       7 entries [r p a K q alpha d], each 1 (log scale),
%                       0 (range) or NaN (automatic), replace that choice.
%     notIdentifiedRatio 0.5: a parameter whose 95% bootstrap interval covers
%                       at least this share of what the search bounds allow
%                       is reported as not identified. Inf: never.
%     weaklyIdentifiedPII 1: a parameter on the log scale whose interval spans
%                       at least this many decades (1 = a factor of 10) is
%                       reported as weakly identified. Inf: never by this rule.
%     tradeoffCorrelation 0.9: pairs of parameters whose rank correlation over
%                       the bootstrap refits is at least this (in size) are
%                       listed as trading off against each other.
%
%   INTERVALS THAT ALLOW FOR MODEL ERROR (spatialWaveModelErrorIntervals.m;
%   revision SW-8). Used by spatialWaveRollingForecast.m. They change scores,
%   not fits: saved fits stay valid.
%     modelErrorIntervals false: as before. true: the rolling forecast also
%                       reports the scores of intervals widened by the
%                       model's own earlier forecast errors.
%     modelErrorWindow  3: past errors at the steps within this many of step h
%                       are pooled for step h.
%     modelErrorMinErrors 10: fewer past errors than this: not widened.
%     modelErrorMaxScale 20: the largest widening factor.
%
%   DAY OF THE WEEK (daily data only)
%     weekdayAdjustment false: as before. true: day-of-week factors are
%                       estimated from the calibration data
%                       (spatialWaveWeekdayFactors.m) and the model curve is
%                       multiplied by them in the fit, the bootstrap and the
%                       forecast. It is a fit setting: a fit saved with the
%                       other value is out of date. Needs smoothfactor1 = 1.
%                       (Revision SW-5.)
%
%   Revision SW-2 (after subepidemicFitSettings.m of the SubEpiPredict toolbox).
%
%   See also fitSpatialWaveModel, spatialWaveSearchSpace, spatialWaveOdeOptions

S = struct( ...
    'logScaleR',true,'logScaleK',true,'logScaleAlpha',true, ...
    'rLB',[],'rUB',[],'pLB',[],'pUB',[],'aLB',[],'aUB',[],'KLB',[],'KUB',[],'qLB',[],'qUB',[], ...
    'rUpperBoundFactor',10, ...
    'countLikelihoodOnRawData',true, ...
    'optimOptions',struct(), ...
    'fitMaxAttempts',3, ...
    'useParallel',false, ...
    'polish',true,'polishGradientTol',10,'polishBudget',100, ...
    'agreementTol',1e-2, ...
    'continuationStart',true, ...
    'firstThresholdStartFactor',3, ...
    'thresholdLowPoints',6, ...
    'bootstrapStarts',1, ...
    'rngSeed',1, ...
    'baselineWindow',[], ...
    'weekdayAdjustment',false, ...
    'PIIscale',[], ...
    'notIdentifiedRatio',0.5, ...
    'weaklyIdentifiedPII',1, ...
    'modelErrorIntervals',false, ...
    'modelErrorWindow',3, ...
    'modelErrorMinErrors',10, ...
    'modelErrorMaxScale',20, ...
    'tradeoffCorrelation',0.9);

if exist('options_fit_settings_sw','file')==2
    U = options_fit_settings_sw();
    if ~(isstruct(U) && isscalar(U))
        error('SpatialWave:FitSettings', ...
            'options_fit_settings_sw must return one struct of settings.');
    end
    names = fieldnames(U);
    for k = 1:numel(names)
        if ~isfield(S,names{k})
            error('SpatialWave:FitSettings', ...
                'Unknown fit setting ''%s''. See spatialWaveFitSettings.m for the list.',names{k});
        end
        S.(names{k}) = U.(names{k});
    end
end

% ---- validation -------------------------------------------------------------
flags = {'logScaleR','logScaleK','logScaleAlpha','countLikelihoodOnRawData','useParallel','polish','continuationStart'};
for k = 1:numel(flags)
    v = S.(flags{k});
    if ~((islogical(v) || isnumeric(v)) && isscalar(v) && (v==0 || v==1))
        error('SpatialWave:FitSettings','Setting %s must be true or false.',flags{k});
    end
    S.(flags{k}) = logical(v);
end
positive = {'polishGradientTol','polishBudget','agreementTol'};
for k = 1:numel(positive)
    v = S.(positive{k});
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>0)
        error('SpatialWave:FitSettings','Setting %s must be a positive number.',positive{k});
    end
end
v = S.rUpperBoundFactor;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=0)
    error('SpatialWave:FitSettings','Setting rUpperBoundFactor must be a number of at least 0.');
end
counts = {'fitMaxAttempts','bootstrapStarts','thresholdLowPoints'};
minimum = [1 0 0];
for k = 1:numel(counts)
    v = S.(counts{k});
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v==fix(v) && v>=minimum(k))
        error('SpatialWave:FitSettings','Setting %s must be an integer of at least %d.',counts{k},minimum(k));
    end
end
bounds = {'rLB','rUB','pLB','pUB','aLB','aUB','KLB','KUB','qLB','qUB'};
for k = 1:numel(bounds)
    v = S.(bounds{k});
    if ~(isempty(v) || (isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v)))
        error('SpatialWave:FitSettings', ...
            'Setting %s must be empty or one finite number (the parameters are shared by all sub-epidemics).',bounds{k});
    end
end
v = S.rngSeed;
if ~(isempty(v) || (isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=0 && v==fix(v) && v<2^32))
    error('SpatialWave:FitSettings','Setting rngSeed must be empty or a nonnegative integer below 2^32.');
end
v = S.notIdentifiedRatio;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && ~isnan(v) && v>0)
    error('SpatialWave:FitSettings','Setting notIdentifiedRatio must be a positive number (Inf: no parameter is reported as not identified).');
end
v = S.modelErrorIntervals;
if ~((islogical(v) || isnumeric(v)) && isscalar(v) && (v==0 || v==1))
    error('SpatialWave:FitSettings','Setting modelErrorIntervals must be true or false.');
end
S.modelErrorIntervals = logical(v);
v = S.modelErrorWindow;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=0 && v==fix(v))
    error('SpatialWave:FitSettings','Setting modelErrorWindow must be a whole number, 0 or more.');
end
v = S.modelErrorMinErrors;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=1 && v==fix(v))
    error('SpatialWave:FitSettings','Setting modelErrorMinErrors must be a positive whole number.');
end
v = S.modelErrorMaxScale;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && ~isnan(v) && v>=1)
    error('SpatialWave:FitSettings','Setting modelErrorMaxScale must be 1 or more.');
end
v = S.weaklyIdentifiedPII;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && ~isnan(v) && v>0)
    error('SpatialWave:FitSettings','Setting weaklyIdentifiedPII must be a positive number (Inf: not used).');
end
v = S.tradeoffCorrelation;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && v>0 && v<=1)
    error('SpatialWave:FitSettings','Setting tradeoffCorrelation must be above 0 and at most 1.');
end
v = S.PIIscale;
if ~(isempty(v) || (isnumeric(v) && isreal(v) && numel(v)==7 && all(isnan(v(:)) | v(:)==0 | v(:)==1)))
    error('SpatialWave:FitSettings', ...
        'Setting PIIscale must be empty or 7 entries [r p a K q alpha d], each 1 (log scale), 0 (range) or NaN (automatic).');
end
v = S.weekdayAdjustment;
if ~((islogical(v) || isnumeric(v)) && isscalar(v) && (v==0 || v==1))
    error('SpatialWave:FitSettings','Setting weekdayAdjustment must be true or false.');
end
S.weekdayAdjustment = logical(v);
v = S.baselineWindow;
if ~(isempty(v) || (isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=1 && v==fix(v)))
    error('SpatialWave:FitSettings','Setting baselineWindow must be empty or a positive integer.');
end
v = S.firstThresholdStartFactor;
if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=1)
    error('SpatialWave:FitSettings','Setting firstThresholdStartFactor must be a number of at least 1 (1: no extra start points).');
end
if ~(isstruct(S.optimOptions) && isscalar(S.optimOptions))
    error('SpatialWave:FitSettings','Setting optimOptions must be a struct of fmincon options.');
end
end
