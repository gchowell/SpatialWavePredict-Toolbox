function U = options_fit_settings_sw_example()
%OPTIONS_FIT_SETTINGS_SW_EXAMPLE Template for optional fit settings.
%   Copy this file to options_fit_settings_sw.m (and rename the function) to
%   override defaults of spatialWaveFitSettings.m. Return only the fields you
%   want to change; every other setting keeps its default. While no
%   options_fit_settings_sw.m exists, the defaults are used.

U = struct();

% ---- speed ------------------------------------------------------------------
% Run the local solves of the search and the bootstrap refits on a parallel
% pool (Parallel Computing Toolbox):
% U.useParallel = true;
%
% Looser fmincon tolerances (faster, slightly less precise optima):
% U.optimOptions = struct('StepTolerance',1e-8,'FunctionTolerance',1e-8);
%
% One MultiStart pass per candidate, and no polish:
% U.fitMaxAttempts = 1;  U.polish = false;

% ---- bounds (empty or absent = the earlier bound, see spatialWaveSearchSpace.m)
% Size of the first sub-epidemic between 20 and 5 million:
% U.KUB = 5e6;
%
% Fix p at 1 (lower bound = upper bound):
% U.pLB = 1;  U.pUB = 1;
%
% The earlier upper bound of r, 5 times the larger of the first two
% observations, without the floor of 10 times the largest count:
% U.rUpperBoundFactor = 0;

% ---- search -------------------------------------------------------------------
% The regular grid of onset thresholds alone, as before revision SW-2:
% U.thresholdLowPoints = 0;
%
% Fit every threshold from the common start set only:
% U.continuationStart = false;

% ---- behaviour ----------------------------------------------------------------
% Fit Poisson and negative-binomial likelihoods to the moving average of span
% smoothfactor1, as the toolbox did before (default: the raw counts):
% U.countLikelihoodOnRawData = false;
%
% More exploration in each bootstrap refit:
% U.bootstrapStarts = 3;

% ---- daily data with a weekly reporting cycle (revision SW-5) -------------------
% Estimate day-of-week factors from the calibration data and multiply the
% model curve by them in the fit, the bootstrap and the forecast. Needs daily
% data (DT = 1) and smoothfactor1 = 1. Fits saved without it are fitted again:
% U.weekdayAdjustment = true;
%
% ---- forecast scores (revision SW-3) ----------------------------------------------
% Number of last observations the flat baseline averages (default: 7 for
% daily data, 1 otherwise). Not a fit setting:
% U.baselineWindow = 14;

% ---- identifiability report (PII; revision SW-6) -----------------------------------
% Report-only settings: saved fits stay valid. A parameter whose 95% bootstrap
% interval covers at least this share of what the search bounds allow is
% reported as not identified (default 0.5):
% U.notIdentifiedRatio = 0.5;
% On the log scale, an interval spanning at least this many decades is
% reported as weakly identified (default 1, a factor of 10):
% U.weaklyIdentifiedPII = 1;
% Pairs of parameters with a rank correlation of at least this size over the
% refits are listed as trading off (default 0.9):
% U.tradeoffCorrelation = 0.9;
% Headline measure per parameter [r p a K q alpha d]: 1 = PII on the log
% scale, 0 = share of the search range, NaN = automatic:
% U.PIIscale = [NaN NaN NaN NaN NaN NaN NaN];

% ---- forecast intervals that allow for model error (revision SW-8) -------------------
% In a rolling forecast, also report the scores of intervals widened by the
% model's own earlier forecast errors. Changes scores only, never a fit:
% U.modelErrorIntervals = true;
% Fewer past errors than this for a step: that step is not widened (default 10):
% U.modelErrorMinErrors = 10;

% ---- repeatable runs ------------------------------------------------------------
% The random number generator is reset with this seed at the start of the fit
% and of the forecast script (default 1). Empty leaves the generator alone:
% U.rngSeed = [];
end
