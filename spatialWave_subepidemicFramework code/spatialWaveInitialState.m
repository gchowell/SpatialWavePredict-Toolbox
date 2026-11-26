function [IC,active0]=spatialWaveInitialState(I0,npatches,onset_fixed)
%SPATIALWAVEINITIALSTATE Initial state and activation mask, one rule for every script.
%   [IC,active0] = spatialWaveInitialState(I0,npatches,onset_fixed)
%
%   onset_fixed = 0 (sub-epidemics start one after another)
%       IC = [I0; 1; ...; 1]. Only the first sub-epidemic is active. The others
%       hold a seed of one case that does not grow until they are activated;
%       spatialWaveIncidence.m removes those seeds from the first incidence.
%   onset_fixed = 1 (all start at the initial time)
%       IC = I0/npatches for every sub-epidemic, all active, no seeds.
%
%   Until revision SW-1 three rules were in use for onset_fixed = 1: the
%   objective optimized a curve from [I0; 1; ...; 1] without removing the
%   seeds, the accepted fit was drawn from I0/npatches, and the forecast used
%   [I0; 1; ...; 1] and removed the seeds. The curve that was scored and
%   bootstrapped was therefore not the curve that had been optimized.
%
%   See also simulateSpatialWave, spatialWaveIncidence

if ~(isnumeric(I0) && isreal(I0) && isscalar(I0) && isfinite(I0) && I0>=0)
    error('SpatialWave:InvalidI0','I0 must be a finite, real, nonnegative scalar.');
end
if ~(isnumeric(npatches) && isreal(npatches) && isscalar(npatches) && ...
        isfinite(npatches) && npatches>=1 && npatches==fix(npatches))
    error('SpatialWave:InvalidPatchCount','npatches must be a positive integer.');
end
if ~(isscalar(onset_fixed) && (onset_fixed==0 || onset_fixed==1))
    error('SpatialWave:InvalidOnsetMode','onset_fixed must be 0 or 1.');
end

IC=ones(npatches,1);
active0=false(npatches,1);
if onset_fixed==0
    IC(1)=I0;
    active0(1)=true;
else
    IC(:)=I0/npatches;
    active0(:)=true;
end
end
