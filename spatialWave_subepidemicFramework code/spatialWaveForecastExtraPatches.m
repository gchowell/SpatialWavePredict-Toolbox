function n=spatialWaveForecastExtraPatches()
%SPATIALWAVEFORECASTEXTRAPATCHES Sub-epidemics a forecast may add beyond the fitted ones.
%   n = spatialWaveForecastExtraPatches() returns 10, the number the forecast
%   script has always added to let the wave travel on ("npatches+10"). It is
%   in one place so that the forecast and the reported size of the wave
%   (plotFit_SW_subepidemicFramework.m) use the same limit. The extra
%   sub-epidemics can start only after the calibration period; see the
%   extension argument of simulateSpatialWave.m.
n=10;
end
