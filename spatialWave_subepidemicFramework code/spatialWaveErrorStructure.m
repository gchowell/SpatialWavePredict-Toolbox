function dist1=spatialWaveErrorStructure(method1,dist1)
%SPATIALWAVEERRORSTRUCTURE The error structure that goes with an estimation method.
%   dist1 = spatialWaveErrorStructure(method1,dist1) returns the value of
%   dist1 every script must use, and stops when the pair cannot run:
%
%     method1 = 0 (least squares)            dist1 = 0, 1 or 2, as given
%     method1 = 2 (Pearson chi-squared)      dist1 = 0, 1 or 2, as given
%     method1 = 1 (Poisson likelihood)       dist1 = 1
%     method1 = 3, 4, 5 (negative binomial)  dist1 = method1
%     method1 = 6 (absolute deviations)      dist1 = 6
%
%   WHY. For methods 3, 4 and 5 the search changed dist1 to the value of
%   method1 while it ran, so the fit files were named '...-method-3-dist-3'.
%   The plotting and forecast scripts read dist1 from options.m; unless the
%   user had also set it there, they looked for '...-method-3-dist-0' and
%   stopped with a missing file. The fit driver and the three scripts now
%   call this function right after options.m, so they agree whatever
%   options.m says about dist1 for these methods.
%
%   See also options, Run_SW_subepidemicFramework

if ~(isnumeric(method1) && isscalar(method1) && any(method1==[0 1 2 3 4 5 6]))
    error('SpatialWave:InvalidSettings','method1 must be 0, 1, 2, 3, 4, 5 or 6.');
end
switch method1
    case {0,2}
        if ~(isnumeric(dist1) && isscalar(dist1) && any(dist1==[0 1 2]))
            error('SpatialWave:InvalidSettings', ...
                ['method1 = %d needs dist1 = 0 (normal), 1 (Poisson) or 2 (negative binomial with ', ...
                 'VAR = factor*mean). Set dist1 in options.m.'],method1);
        end
    case 1
        dist1=1;
    case {3,4,5}
        dist1=method1;
    case 6
        dist1=6;
end
end
