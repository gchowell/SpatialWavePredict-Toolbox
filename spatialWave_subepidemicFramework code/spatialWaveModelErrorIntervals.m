function [adjusted,scales] = spatialWaveModelErrorIntervals(forecasts,DT,settings)
%SPATIALWAVEMODELERRORINTERVALS Widen forecast intervals by the model's own past errors.
%   [adjusted,scales] = spatialWaveModelErrorIntervals(forecasts,DT,settings)
%
%   WHY. The prediction intervals of the toolbox come from the bootstrap
%   refits and the error structure. Both take the model as right. When the
%   epidemic does something the model cannot represent, the interval stays
%   narrow around a forecast that is in the wrong place. This function
%   measures how far the model's earlier forecasts were from what was then
%   observed and widens the intervals of a new forecast accordingly. The
%   median forecast is not changed.
%
%   forecasts   one row per forecast date, model and step:
%               [date number, model, step h, observation, 23 quantiles]
%               (R.forecasts of spatialWaveRollingForecast.m; quantile
%               levels as in spatialWaveQuantileScores.m)
%   DT          days per time unit of the data (1, 7 or 365)
%   settings    spatialWaveFitSettings.m ([] = the settings in force). Used:
%               modelErrorWindow, modelErrorMinErrors, modelErrorMaxScale.
%
%   adjusted    the same rows with the quantiles widened; rows that could
%               not be adjusted are returned as they were
%   scales      one row per row of forecasts: [date number, model, step h,
%               scale, number of past errors used]. Scale 1 = not widened.
%
%   HOW, for the forecast of model m from date D at step h:
%   1. Past errors. Take the forecasts of the SAME model from EARLIER dates
%      d < D, at the steps j within modelErrorWindow (3) of h, and only where
%      the observation was already known on date D (d + j time units <= D).
%      Nothing from date D on is used.
%   2. Each past error is measured against the half-width of its own 95%
%      interval: (y - median)/(upper - median) for an observation above the
%      median, (median - y)/(median - lower) below it. A value of 1 means the
%      observation sat exactly on the interval limit.
%   3. The scale is the 95% quantile of those values, at least 1 (intervals
%      are never narrowed) and at most modelErrorMaxScale (20). With fewer
%      than modelErrorMinErrors (10) past errors the scale is 1.
%   4. Every quantile is moved away from the median by the scale:
%         q' = median + scale*(q - median), not below 0.
%
%   A model whose intervals were right in the past gets a scale of 1: about
%   95% of its past errors are below 1.
%
%   WHAT IT DOES NOT DO. It does not move the forecast. A forecast in the
%   wrong place stays there; the interval then says so. It needs earlier
%   forecast dates, so the first dates of a run are not adjusted. Past errors
%   of overlapping forecasts are not independent, and with few dates the
%   scale is rough. The flat baseline (model 0) is left as it is: its
%   intervals are already built from past errors.
%
%   Revision SW-8.
%
%   See also spatialWaveQuantileScores, spatialWaveRollingForecast

if nargin<3 || isempty(settings)
    settings=spatialWaveFitSettings();
end
window=localField(settings,'modelErrorWindow',3);
minErrors=localField(settings,'modelErrorMinErrors',10);
maxScale=localField(settings,'modelErrorMaxScale',20);
if size(forecasts,2)~=27
    error('SpatialWave:ModelErrorIntervals', ...
        'forecasts must have 27 columns: date number, model, step, observation and 23 quantiles.');
end
if ~(isnumeric(DT) && isscalar(DT) && DT>0)
    error('SpatialWave:ModelErrorIntervals','DT must be the number of days per time unit.');
end
n=size(forecasts,1);
adjusted=forecasts;
scales=[forecasts(:,1:3) ones(n,1) zeros(n,1)];

% normalised error of every row (step 2)
median1=forecasts(:,4+12);
lower1=forecasts(:,4+2);
upper1=forecasts(:,4+22);
y=forecasts(:,4);
normalised=NaN(n,1);
for i=1:n
    if ~(isfinite(y(i)) && isfinite(median1(i)) && isfinite(lower1(i)) && isfinite(upper1(i)))
        continue
    end
    if y(i)>=median1(i)
        distance=y(i)-median1(i); halfWidth=upper1(i)-median1(i);
    else
        distance=median1(i)-y(i); halfWidth=median1(i)-lower1(i);
    end
    if distance==0
        normalised(i)=0;
    elseif halfWidth>0
        normalised(i)=distance/halfWidth;
    else
        normalised(i)=Inf;      % an interval of no width that missed
    end
end
knownOn=forecasts(:,1)+forecasts(:,3)*DT;   % the day the observation of a row became known

models=unique(forecasts(:,2))';
for m=models(models~=0)
    rowsModel=find(forecasts(:,2)==m);
    for i=rowsModel'
        D=forecasts(i,1);
        h=forecasts(i,3);
        past=rowsModel(forecasts(rowsModel,1)<D & knownOn(rowsModel)<=D & ...
            abs(forecasts(rowsModel,3)-h)<=window & ~isnan(normalised(rowsModel)));
        scales(i,5)=numel(past);
        if numel(past)<minErrors || ~all(isfinite(forecasts(i,5:27)))
            continue
        end
        scale=min(max(quantile(normalised(past),0.95),1),maxScale);
        scales(i,4)=scale;
        adjusted(i,5:27)=max(median1(i)+scale*(forecasts(i,5:27)-median1(i)),0);
    end
end
end

% -----------------------------------------------------------------------------
function v = localField(S,name,default)
if isstruct(S) && isfield(S,name)
    v = S.(name);
else
    v = default;
end
end
