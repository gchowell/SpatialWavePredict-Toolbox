function curve = spatialWaveWeekdayCurve(weekday,n)
%SPATIALWAVEWEEKDAYCURVE Day-of-week factor of each of the first n time points.
%   curve = spatialWaveWeekdayCurve(weekday,n) returns an n-by-1 vector:
%   element t is the factor of the day of the week of time point t, with
%   time point 1 the first observation of the calibration series the factors
%   were estimated from (spatialWaveWeekdayFactors.m). The 7 factors repeat,
%   so n may reach beyond the calibration period (forecasts).
%
%   Returns [] when no adjustment applies. Callers leave the model curve as
%   it is in that case, so fits without the adjustment are computed exactly
%   as before revision SW-5.
%
%   Revision SW-5.
if ~(isstruct(weekday) && isfield(weekday,'applied') && weekday.applied)
    curve = [];
    return
end
if ~(isnumeric(n) && isscalar(n) && n>=0 && n==fix(n))
    error('SpatialWave:Weekday','n must be a whole number, 0 or more.');
end
curve = weekday.byPosition(mod((0:n-1)',7)+1);
end
