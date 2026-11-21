function objfunction=plotModifiedLogisticGrowthPatchMethods1(z)
%PLOTMODIFIEDLOGISTICGROWTHPATCHMETHODS1 Fitting objective of the spatial wave model.
%   objfunction = plotModifiedLogisticGrowthPatchMethods1(z), with
%   z = [r p a K q alpha d]. The model settings and the data are read from
%   the global variables the fitting drivers set (flag1, method1, timevect,
%   ydata, I0, npatches, onset_thr, typedecline1, onset_fixed).
%
%   Revision SW-1. The curve comes from simulateSpatialWave.m, the simulator
%   every other script uses, so the curve that is optimized is the curve
%   that is scored, bootstrapped, plotted and forecast. A parameter vector
%   whose curve cannot be simulated (the solver stops early, or the state
%   is not finite or not real) gets the penalty 1e10 the objective already
%   used for an all-zero curve; before, such a vector could stop the fit
%   with an error or return a complex value.
%
%   See also simulateSpatialWave, getSpatialWaveObjective

global flag1 method1 timevect ydata I0 npatches onset_thr typedecline1 yfit

global onset_fixed

r=z(1);
p=z(2);
a=z(3);
K=z(4);
q=z(5);

alpha=z(end-1);
d=z(end);

try
    [~,~,totinc]=simulateSpatialWave(timevect,I0,npatches,onset_fixed,onset_thr,flag1,typedecline1,r,p,a,K,q);
catch simulationError
    if strncmp(simulationError.identifier,'SpatialWave:',12) && ...
            any(strcmp(simulationError.identifier,{'SpatialWave:IncompleteIntegration', ...
            'SpatialWave:InvalidSegment','SpatialWave:InvalidTrajectory', ...
            'SpatialWave:InvalidGrowthParameters','SpatialWave:ActivationProgressFailure'}))
        objfunction=10^10;
        return
    end
    rethrow(simulationError);
end

yfit=totinc;

[objfunction,yfit]=getSpatialWaveObjective(ydata,yfit,method1,alpha,d);

if ~(isscalar(objfunction) && isreal(objfunction) && isfinite(objfunction))
    objfunction=10^10;
end
end
