function scores = spatialWaveForecastScores(forecasts,original)
%SPATIALWAVEFORECASTSCORES Scores over the first h steps from quantile forecasts.
%   scores = spatialWaveForecastScores(forecasts,original)
%
%   forecasts   rows [date number, model, step h, observation, 23 quantiles]
%               (R.forecasts of spatialWaveRollingForecast.m, or the first
%               output of spatialWaveModelErrorIntervals.m)
%   original    the scores of the same forecasts as the forecast script gave
%               them (R.scores): rows [date number, model, h, MAE, MSE,
%               coverage, WIS, WIS of the baseline]
%
%   scores      the layout of original. MAE and MSE are copied from original
%               (they belong to the point forecast, which the quantiles do
%               not change); coverage and WIS are computed from the
%               quantiles (spatialWaveQuantileScores.m) and averaged over the
%               first h steps; the last column is the WIS of the baseline
%               (model 0) of the same date and horizon.
%
%   Called with R.forecasts it reproduces R.scores. Called with widened
%   forecasts it gives the scores of the widened intervals.
%
%   Revision SW-8.
scores=original;
scores(:,6:8)=NaN;
keys=unique(forecasts(:,1:2),'rows');
for k=1:size(keys,1)
    rows=find(forecasts(:,1)==keys(k,1) & forecasts(:,2)==keys(k,2));
    [~,order]=sort(forecasts(rows,3));
    rows=rows(order);
    steps=forecasts(rows,3);
    if ~isequal(steps(:)',1:numel(steps))
        error('SpatialWave:ForecastScores','The steps of date %g, model %g are not 1, 2, ..., H.',keys(k,1),keys(k,2));
    end
    [wisStep,insideStep]=spatialWaveQuantileScores(forecasts(rows,5:27),forecasts(rows,4));
    wisOver=cumsum(wisStep)./(1:numel(steps))';
    coverOver=100*cumsum(insideStep)./(1:numel(steps))';
    for j=1:numel(steps)
        target=scores(:,1)==keys(k,1) & scores(:,2)==keys(k,2) & scores(:,3)==steps(j);
        scores(target,6)=coverOver(j);
        scores(target,7)=wisOver(j);
    end
end
% the WIS of the baseline next to every row
for i=1:size(scores,1)
    base=scores(:,1)==scores(i,1) & scores(:,2)==0 & scores(:,3)==scores(i,3);
    if any(base)
        scores(i,8)=scores(find(base,1),7);
    end
end
end
