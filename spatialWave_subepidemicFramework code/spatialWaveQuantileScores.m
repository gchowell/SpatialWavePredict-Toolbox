function [wisStep,insideStep,levels] = spatialWaveQuantileScores(quantiles,observed)
%SPATIALWAVEQUANTILESCORES WIS and 95% coverage of quantile forecasts, step by step.
%   [wisStep,insideStep,levels] = spatialWaveQuantileScores(quantiles,observed)
%
%   quantiles   H-by-23 forecast quantiles, one row per forecast step, at the
%               levels 0.01, 0.025, 0.05, 0.10, ..., 0.95, 0.975, 0.99 (the
%               fourth output of computeWIS.m; returned here as levels)
%   observed    H observations
%
%   wisStep     weighted interval score of each step, with the intervals and
%               weights of computeWIS.m: 11 central intervals (alpha = 0.02,
%               0.05, 0.1, ..., 0.9) and the median.
%   insideStep  1 where the observation lies in the 95% interval (0.025 to
%               0.975), 0 where it does not.
%
%   The scores computeWIS.m and computeforecastperformance.m report for a
%   horizon h are averages over the first h steps:
%       WIS over h      = mean(wisStep(1:h))
%       coverage over h = 100*mean(insideStep(1:h))
%   Computed from the quantiles those functions use, the two agree with them.
%   A step with a quantile that is not finite gets NaN.
%
%   Revision SW-8. Used by spatialWaveModelErrorIntervals.m to score
%   forecasts whose intervals have been widened.

levels=[0.010 0.025 0.050 0.100 0.150 0.200 0.250 0.300 0.350 0.400 0.450 0.500 ...
    0.550 0.600 0.650 0.700 0.750 0.800 0.850 0.900 0.950 0.975 0.990];
if size(quantiles,2)~=numel(levels)
    error('SpatialWave:QuantileScores','quantiles must have 23 columns, one per level (see the header).');
end
observed=observed(:);
H=size(quantiles,1);
if numel(observed)~=H
    error('SpatialWave:QuantileScores','One observation is needed for every row of quantiles.');
end
alphas=[0.02 0.05 0.1:0.1:0.9];
K=numel(alphas);
wisStep=NaN(H,1);
insideStep=NaN(H,1);
for i=1:H
    q=quantiles(i,:);
    y=observed(i);
    if ~(all(isfinite(q)) && isfinite(y))
        continue
    end
    total=0.5*abs(y-q(12));
    for k=1:K
        lower=q(k);            % level alpha/2
        upper=q(24-k);         % level 1-alpha/2
        score=(upper-lower)+(2/alphas(k))*(lower-y)*(y<lower)+(2/alphas(k))*(y-upper)*(y>upper);
        total=total+(alphas(k)/2)*score;
    end
    wisStep(i)=total/(K+0.5);
    insideStep(i)=double(y>=q(2) && y<=q(22));
end
end
