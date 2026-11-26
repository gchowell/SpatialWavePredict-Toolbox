function [value,detail] = spatialWaveObjective(z,context)
%SPATIALWAVEOBJECTIVE Fitting objective with an explicit context, no globals.
%   [value,detail] = spatialWaveObjective(z,context)
%
%   z        parameter vector [r p a K q alpha d]
%   context  struct with fields
%              timevect, ydata   observation times and the series fitted
%              I0                initial value (first observation)
%              npatches          number of sub-epidemics the model may have
%              onset_fixed       0 one after another, 1 all at the first time
%              onset_thr         onset threshold C_thr
%              flag1             growth model, 0 to 5
%              typedecline1      1 exponential, 2 power-law decline of sizes
%              method1           0 least squares, 1 Poisson, 2 Pearson
%                                chi-squared, 3/4/5 negative binomial,
%                                6 sum of absolute deviations
%              odeOptions        odeset structure (spatialWaveOdeOptions)
%
%   value    the number the earlier objective
%            (plotModifiedLogisticGrowthPatchMethods1.m) returns for the same
%            model: SSE, the negative log-likelihood without its data-only
%            constant, or the sum of absolute deviations. AICc is computed
%            from it as before (getSpatialWaveObjective.m does the arithmetic
%            for both).
%   detail   struct: success, message, nActive (sub-epidemics active by the
%            end of the time grid), activation (fifth output of
%            simulateSpatialWave), curve (the fitted incidence).
%
%   A parameter vector is INVALID (detail.success = false, value = 1e10)
%   when the simulation fails, the curve is identically zero, or the
%   objective is not a finite real number.
%
%   Because everything the objective needs is in CONTEXT, a handle to it can
%   be sent to the workers of a parallel pool. The earlier objective read
%   global variables, which workers do not share.
%
%   Revision SW-2 (after subepidemicObjective.m of the SubEpiPredict toolbox).
%
%   See also fitSpatialWaveModel, simulateSpatialWave, getSpatialWaveObjective

value = 1e10;
detail = struct('success',false,'message','','nActive',0,'activation',[],'curve',[]);
% A search visits parameter vectors whose curve cannot be integrated; they
% get the penalty below. ode15s would print a warning for each of them, on
% the client or on a worker, so that warning is switched off for this call.
odeWarning = warning('off','MATLAB:ode15s:IntegrationTolNotMet');
try
    if ~(isnumeric(z) && isreal(z) && numel(z)==7 && all(isfinite(z(:))))
        detail.message = 'Parameter vector must hold 7 finite real values [r p a K q alpha d].';
        warning(odeWarning);
        return
    end
    [~,~,yfit,~,activation] = simulateSpatialWave(context.timevect,context.I0,context.npatches, ...
        context.onset_fixed,context.onset_thr,context.flag1,context.typedecline1, ...
        z(1),z(2),z(3),z(4),z(5),context.odeOptions);
    detail.activation = activation;
    detail.nActive = activation.nActive;
    % Revision SW-5: day-of-week factors (spatialWaveWeekdayCurve.m). The
    % field is absent or empty when no adjustment applies.
    if isfield(context,'weekdayFactor') && ~isempty(context.weekdayFactor)
        yfit = yfit(:).*context.weekdayFactor(:);
    end
    if sum(yfit)==0
        detail.message = 'The fitted curve is identically zero.';
        warning(odeWarning);
        return
    end
    [v,curve] = getSpatialWaveObjective(context.ydata(:),yfit(:),context.method1,z(6),z(7));
    if ~(isscalar(v) && isreal(v) && isfinite(v))
        detail.message = 'The objective is not a finite real number at this parameter vector.';
        warning(odeWarning);
        return
    end
    value = v;
    detail.success = true;
    detail.curve = curve;
catch ME
    value = 1e10;
    detail.success = false;
    detail.message = sprintf('%s: %s',ME.identifier,ME.message);
end
warning(odeWarning);
end
