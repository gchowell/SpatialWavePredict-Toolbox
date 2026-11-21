function [totinc,patchIncidence]=spatialWaveIncidence(x,onset_fixed)
%SPATIALWAVEINCIDENCE Incidence per time interval from the cumulative states.
%   [totinc,patchIncidence] = spatialWaveIncidence(x,onset_fixed)
%
%   x holds the states returned by simulateSpatialWave.m: rows are the
%   observation times, starting at the initial time; columns are the
%   sub-epidemics. Row 1 of the incidence is the initial state and each later
%   row is the increase over one time interval. For onset_fixed = 0 the seed
%   of one case of sub-epidemics 2..n is removed from row 1, so totinc(1) is
%   the first observation I0. sum(patchIncidence,2) equals totinc.
%
%   See also simulateSpatialWave, spatialWaveInitialState

if ~(isnumeric(x) && ismatrix(x) && ~isempty(x))
    error('SpatialWave:InvalidTrajectory','x must be a nonempty numeric matrix.');
end
if ~(isscalar(onset_fixed) && (onset_fixed==0 || onset_fixed==1))
    error('SpatialWave:InvalidOnsetMode','onset_fixed must be 0 or 1.');
end
patchIncidence=[x(1,:);diff(x,1,1)];
if onset_fixed==0
    patchIncidence(1,2:end)=patchIncidence(1,2:end)-1;
end
totinc=sum(patchIncidence,2);
end
