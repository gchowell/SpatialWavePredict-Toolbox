function pairs = spatialWaveSameCurveModels(fits,tolerance)
%SPATIALWAVESAMECURVEMODELS Ranked models with the same curve over the calibration period.
%   pairs = spatialWaveSameCurveModels(fits) takes the cell array of accepted
%   fits of a final ranking (finalRanking.fits) and returns one row [i j],
%   i < j, for every two ranked models whose fitted curves agree to within
%   TOLERANCE (default 1e-6) of the largest fitted value.
%
%   WHY. A wave candidate whose first sub-epidemic has not reached the onset
%   threshold by the end of the calibration period has ONE active
%   sub-epidemic. Its curve, likelihood and AICc do not depend on the
%   threshold or on the decline rate q, so the same fit appears once for
%   every threshold that was not reached. Such models are different only in
%   the forecast, through a threshold and a q that the data do not inform.
%   The fit driver reports them (warning SpatialWave:SameCurveModels); it
%   does not remove them.
%
%   Revision SW-2.
%
%   See also Run_SW_subepidemicFramework, finalizeSpatialWaveFits

if nargin<2 || isempty(tolerance)
    tolerance = 1e-6;
end
pairs = zeros(0,2);
n = numel(fits);
for i = 1:n-1
    a = fits{i}.bestfit(:);
    for j = i+1:n
        b = fits{j}.bestfit(:);
        if numel(a)==numel(b) && max(abs(a-b))<=tolerance*max(1,max(abs(a)))
            pairs(end+1,:) = [i j]; %#ok<AGROW>
        end
    end
end
end
