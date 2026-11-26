function options=spatialWaveOdeOptions()
%SPATIALWAVEODEOPTIONS One place for the ODE solver tolerances of the toolbox.
%   options = spatialWaveOdeOptions() returns the odeset structure that
%   simulateSpatialWave.m uses whenever no solver options are passed, so the
%   fit, the bootstrap, the plots and the forecast share one setting.
%
%   RelTol 1e-8 and AbsTol 1e-10, the values of subepidemicOdeOptions.m in the
%   SubEpiPredict toolbox. Incidence is a difference of cumulative values that
%   are much larger than it, so the odeset defaults (1e-3, 1e-6) leave an
%   error in the fitted incidence that is visible in the AICc and makes the
%   finite-difference gradients of fmincon noisy.
%
%   Do not set Events here: simulateSpatialWave.m owns them.
options=odeset('RelTol',1e-8,'AbsTol',1e-10);
end
