function patchCap=spatialWavePatchCap(nActive,onset_thr,onset_fixed,npatches_fixed)
%SPATIALWAVEPATCHCAP Number of sub-epidemics a ranked candidate was fitted with.
%   patchCap = spatialWavePatchCap(nActive,onset_thr,onset_fixed,npatches_fixed)
%
%   The search fits every candidate with a CONFIGURED number of sub-epidemics
%   (1 or npatches_fixed when onset_fixed = 0; 1 to npatches_fixed when
%   onset_fixed = 1) and stores in column 1 of its ranking the number that
%   became ACTIVE during the calibration period, which can be smaller. This
%   function recovers the configured number from a row of the ranking:
%
%     onset_fixed = 1   all sub-epidemics start at once: configured = active
%     onset_fixed = 0   one active sub-epidemic with threshold 0 is the
%                       single sub-epidemic candidate (a wave with threshold
%                       0 starts all its sub-epidemics at once, so it has
%                       more than one active); every other row is a wave
%                       fitted with npatches_fixed
%
%   WHY IT MATTERS. Until revision SW-1 the selected-model refit, the
%   bootstrap and the forecast used the ACTIVE count as the size of the
%   model. A wave candidate with one active sub-epidemic was then refitted
%   as a single sub-epidemic model: no bootstrap dataset and no forecast
%   could start the second sub-epidemic that the fitted threshold allows.
%   The active count still enters the AICc, as before.
%
%   See also fittingModifiedLogisticFunctionPatchABC,
%            fittingModifiedLogisticFunctionPatchMultiple, simulateSpatialWave

if ~(isscalar(nActive) && isfinite(nActive) && nActive>=1 && nActive==fix(nActive))
    error('SpatialWave:InvalidPatchCount','The active count must be a positive integer.');
end
if ~(isscalar(npatches_fixed) && isfinite(npatches_fixed) && npatches_fixed>=1 && npatches_fixed==fix(npatches_fixed))
    error('SpatialWave:InvalidPatchCount','npatches_fixed must be a positive integer.');
end
if onset_fixed==1 || (nActive==1 && onset_thr==0)
    patchCap=nActive;
else
    patchCap=max(npatches_fixed,nActive);
end
patchCap=full(double(patchCap));
end
