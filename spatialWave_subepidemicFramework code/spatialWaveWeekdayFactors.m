function weekday = spatialWaveWeekdayFactors(y,settings,seriesInfo)
%SPATIALWAVEWEEKDAYFACTORS Day-of-week factors of a daily series.
%   weekday = spatialWaveWeekdayFactors(y,settings)
%   weekday = spatialWaveWeekdayFactors(y,settings,seriesInfo)
%
%   y           the calibration series (daily incidence), a vector.
%   settings    the fit settings (spatialWaveFitSettings.m). The factors are
%               estimated only when settings.weekdayAdjustment is true; a
%               struct without that field (a fit saved before revision SW-5)
%               counts as false.
%   seriesInfo  optional second output of spatialWaveFitSeries.m. The
%               adjustment needs the model to be fitted to the series as
%               reported; with a moving average (smoothfactor1 above 1) the
%               function stops.
%
%   weekday.applied      true when factors were estimated.
%   weekday.byPosition   7 factors with mean 1. Factor k belongs to the
%                        observations 1+(k-1), 8+(k-1), 15+(k-1), ... of y:
%                        position 1 is the day of the week of the first
%                        observation. All 1 when not applied.
%   weekday.nRatios      number of days each factor rests on.
%   weekday.reason       why the factors were not estimated ('' when applied).
%
%   HOW. Each observation is divided by the mean of the 7 days centred on it
%   (the first and last 3 days have no such mean). The factor of a day of the
%   week is the MEDIAN of its ratios, so a single holiday does not move it.
%   The 7 factors are divided by their mean. A factor below 0.05 is set to
%   0.05 before that (a day of the week without reports).
%
%   USE. The model curve is multiplied by the factor of each day
%   (spatialWaveWeekdayCurve.m) wherever it is compared with data or turned
%   into a forecast: in the objective, the ranked fit, the bootstrap, the
%   fit plot and the forecast. The sub-epidemic model itself describes the
%   series without the weekly cycle.
%
%   The factors are fixed numbers taken from the data before the fit. They
%   are not counted as parameters in the AICc, so the AICc of a fit with the
%   adjustment must not be compared with that of a fit without it.
%
%   Needs at least 21 observations and 2 ratios for every day of the week;
%   otherwise no adjustment is made and a warning says so.
%
%   Revision SW-5.

weekday = struct('applied',false,'byPosition',ones(7,1),'nRatios',zeros(7,1),'reason','');
if nargin<2 || ~isstruct(settings) || ~isfield(settings,'weekdayAdjustment') || ~isequal(logical(settings.weekdayAdjustment),true)
    weekday.reason = 'weekdayAdjustment is not switched on';
    return
end
if ~(isnumeric(y) && isreal(y) && isvector(y) && all(isfinite(y(:))))
    error('SpatialWave:Weekday','The series must be a finite real vector.');
end
if nargin>=3 && isstruct(seriesInfo) && isfield(seriesInfo,'usedRawData') && ~seriesInfo.usedRawData
    error('SpatialWave:WeekdayWithSmoothing', ...
        ['weekdayAdjustment = true cannot be combined with a moving average of the calibration data ', ...
         '(smoothfactor1 above 1): the moving average has already removed the weekly cycle. ', ...
         'Set smoothfactor1 = 1 in options.m or weekdayAdjustment = false in options_fit_settings_sw.m.']);
end
y = double(y(:));
n = numel(y);
if n<21
    weekday.reason = sprintf('the series has %d observations; 21 are needed',n);
    warning('SpatialWave:WeekdayNotEstimated','No day-of-week adjustment: %s.',weekday.reason);
    return
end
ratio = NaN(n,1);
for t = 4:n-3
    centre = mean(y(t-3:t+3));
    if centre>0
        ratio(t) = y(t)/centre;
    end
end
position = mod((0:n-1)',7)+1;
factors = ones(7,1);
for k = 1:7
    mine = ratio(position==k & isfinite(ratio));
    weekday.nRatios(k) = numel(mine);
    if numel(mine)>=2
        factors(k) = median(mine);
    end
end
if any(weekday.nRatios<2)
    weekday.reason = 'fewer than 2 usable weeks for a day of the week (too many zero counts)';
    warning('SpatialWave:WeekdayNotEstimated','No day-of-week adjustment: %s.',weekday.reason);
    return
end
factors = max(factors,0.05);
weekday.byPosition = factors/mean(factors);
weekday.applied = true;
end
