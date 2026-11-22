function tag = spatialWaveOutputTag(nSubepidemics,smoothing,calibrationLength)
%SPATIALWAVEOUTPUTTAG The part of an output file name that identifies the fit.
%   tag = spatialWaveOutputTag() returns, for the settings in force,
%       'npatchesfixed-<n>-smoothing-<s>-calibrationperiod-<c>-'
%   from npatches_fixed, smoothfactor1 and calibrationperiod1 (the global
%   variables options.m sets).
%   tag = spatialWaveOutputTag(nSubepidemics,smoothing,calibrationLength)
%   builds it from the three values given.
%
%   WHY. The files of a fit (modifiedLogisticPatch-..., FinalFit-..., ABC-...)
%   have always carried these three settings in their names. The tables
%   written from a fit (ranked(k)-..., quantile-..., performance-...,
%   parameters-..., identifiability-..., Rt-..., the rolling forecast scores)
%   did not, or carried only some. Two runs that differed only in the
%   calibration window, the smoothing or the number of sub-epidemics
%   therefore wrote the same file names, and the later run replaced the
%   tables of the earlier one. With this tag in every name they no longer do.
%
%   Tables written before this revision keep their old names in ./output and
%   are not touched.
%
%   Revision SW-9.
if nargin<3
    global npatches_fixed smoothfactor1 calibrationperiod1 %#ok<GVMIS,TLEV>
    nSubepidemics=npatches_fixed;
    smoothing=smoothfactor1;
    calibrationLength=calibrationperiod1;
end
values={nSubepidemics,smoothing,calibrationLength};
names={'npatches_fixed','smoothfactor1','calibrationperiod1'};
for k=1:3
    if ~(isnumeric(values{k}) && isscalar(values{k}) && isfinite(values{k}))
        error('SpatialWave:OutputTag', ...
            '%s is not set. Call options.m (or a toolbox entry point) before writing output files.',names{k});
    end
end
tag=sprintf('npatchesfixed-%s-smoothing-%s-calibrationperiod-%s-',num2str(nSubepidemics),num2str(smoothing),num2str(calibrationLength));
end
