function spatialWaveSkippedRefits(label,rank1,nSkipped,nRefits,firstMessage)
%SPATIALWAVESKIPPEDREFITS One policy for bootstrap refits that cannot be simulated.
%   spatialWaveSkippedRefits(label,rank1,nSkipped,nRefits,firstMessage)
%
%   A forecast extends every bootstrap refit past the calibration period; an
%   extreme parameter set can fail there. nSkipped of the nRefits of ranked
%   model rank1 were left out. Then (the policy of subepidemicSkippedRefits.m
%   in the SubEpiPredict toolbox, with its default limit of 5%):
%     none skipped       nothing happens
%     all skipped        error SpatialWave:NoForecastCurves
%     more than 5%       error SpatialWave:TooManyRefitsSkipped
%     otherwise          warning SpatialWave:ForecastCurvesSkipped
%   The refits that fail tend to be the extreme ones, so leaving many out
%   would narrow the interval; that is why there is a limit. Before, such
%   refits were dropped without any message.

maxFraction=0.05;
if nSkipped==0
    return
end
if nargin<5 || isempty(firstMessage)
    firstMessage='(no message)';
end
if nSkipped>=nRefits
    error('SpatialWave:NoForecastCurves', ...
        '%s, ranked model %d: none of the %d bootstrap refits could be simulated. First failure: %s', ...
        label,rank1,nRefits,firstMessage);
elseif nSkipped/nRefits>maxFraction
    error('SpatialWave:TooManyRefitsSkipped', ...
        ['%s, ranked model %d: %d of %d bootstrap refits could not be simulated, more than %g%%. ', ...
         'First failure: %s'],label,rank1,nSkipped,nRefits,100*maxFraction,firstMessage);
else
    warning('SpatialWave:ForecastCurvesSkipped', ...
        '%s, ranked model %d: %d of %d bootstrap refits could not be simulated and were left out. First failure: %s', ...
        label,rank1,nSkipped,nRefits,firstMessage);
end
end
