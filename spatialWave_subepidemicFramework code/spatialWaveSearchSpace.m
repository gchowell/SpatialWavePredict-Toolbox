function [LB,UB,isLog,z0] = spatialWaveSearchSpace(data,flag1,method1,settings)
%SPATIALWAVESEARCHSPACE Bounds, coordinate scales and seed of the fit.
%   [LB,UB,isLog,z0] = spatialWaveSearchSpace(data,flag1,method1,settings)
%
%   Layout of every vector: [r p a K q alpha d], 7 entries.
%
%   data      calibration incidence (the unsmoothed series), a vector
%   flag1     growth model, 0 to 5
%   method1   estimation method, 0 to 6
%   settings  spatialWaveFitSettings.m (optional)
%
%   LB, UB  The bounds the toolbox used before (getSpatialWaveBounds.m, and
%           for alpha and d the bounds of the search), with one exception:
%           the upper bound of r is never below settings.rUpperBoundFactor
%           (default 10) times the largest count. The earlier bound, 5 times
%           the larger of the first two observations, depends on how the
%           series happens to start (5 for a series that starts with 1, 0)
%           and can be reached by a fit with a small p, which needs a large
%           r. rUpperBoundFactor = 0 gives the earlier bounds exactly. The
%           Gompertz bounds of r (1e-4 to 1e3) are not changed.
%           The settings rLB ... qUB replace a bound; a lower bound equal to
%           its upper bound fixes the parameter. Reversed bounds are an error.
%   isLog   true for the coordinates searched on a log scale: r, K and
%           alpha, when the setting is on, the lower bound is positive and
%           the parameter is free. With K between 20 and 1e8 on a linear
%           scale, 20 of 21 random start points lie above 5e6.
%   z0      the seed of the earlier search (r = 0.2, p = 0.9, a = 1, K = sum
%           of the data, q = 0.3, alpha = 1, d = 1), projected into the bounds.
%
%   No global variables are used.
%
%   Revision SW-2 (after subepidemicSearchSpace.m of the SubEpiPredict toolbox).
%
%   See also spatialWaveFitSettings, fitSpatialWaveModel, getSpatialWaveBounds

if nargin<4 || isempty(settings)
    settings = spatialWaveFitSettings();
end
if ~(isnumeric(data) && isreal(data) && isvector(data) && ~isempty(data) && all(isfinite(data(:))))
    error('SpatialWave:InvalidData','data must be a nonempty finite real vector.');
end
if ~(isnumeric(method1) && isscalar(method1) && any(method1==0:6))
    error('SpatialWave:InvalidMethod','method1 must be an integer from 0 to 6.');
end
data = double(data(:));

% bounds of [alpha d] by estimation method
switch method1
    case {0,1,2,6}
        LBe = [0 0];      UBe = [0 0];
    case {3,4}
        LBe = [1e-8 0];   UBe = [1e5 0];
    case 5
        LBe = [1e-8 0.2]; UBe = [1e5 1e2];
end
[LB,UB] = getSpatialWaveBounds(flag1,data,LBe,UBe);

% the upper bound of r follows the size of the counts
if flag1~=5 && settings.rUpperBoundFactor>0
    UB(1) = max(UB(1),settings.rUpperBoundFactor*max(abs(data)));
end

% bounds set by the user
names = {'r','p','a','K','q'};
for k = 1:numel(names)
    lower = settings.([names{k} 'LB']);
    upper = settings.([names{k} 'UB']);
    if ~isempty(lower), LB(k) = lower(1); end
    if ~isempty(upper), UB(k) = upper(1); end
end
if any(~isfinite(LB)) || any(~isfinite(UB)) || any(LB>UB)
    error('SpatialWave:ReversedBounds', ...
        'Bounds must be finite with lower <= upper. Check rLB ... qUB in options_fit_settings_sw.m.');
end

isLog = false(1,7);
isLog(1) = settings.logScaleR;
isLog(4) = settings.logScaleK;
isLog(6) = settings.logScaleAlpha;
isLog = isLog & LB>0 & UB>LB;

z0 = [0.2 0.9 1 sum(data) 0.3 1 1];
z0 = min(max(z0,LB),UB);
end
