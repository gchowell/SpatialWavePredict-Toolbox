function onset_thrs = spatialWaveThresholdGrid(data,smoothfactor1,I0,nLow)
%SPATIALWAVETHRESHOLDGRID Onset thresholds C_thr the search tries.
%   onset_thrs = spatialWaveThresholdGrid(data,smoothfactor1,I0,nLow)
%
%   data           calibration incidence, a vector with n values
%   smoothfactor1  span of the moving average that shapes the grid
%   I0             first observation
%   nLow           thresholds added at the low end (setting
%                  thresholdLowPoints; default 6, 0 for none)
%
%   REGULAR GRID (as before). With C the cumulative sum of the smoothed
%   series: 0, then C(1), then n-2 further values evenly spaced up to (not
%   including) C(n).
%
%   LOW END (revision SW-2). The regular grid is even in cumulative cases, so
%   its first step is about C(n)/n: 18,400 cases for the bundled COVID-19
%   series. A wave whose second sub-epidemic starts earlier than that was not
%   among the candidates, and the best model of that series sat on the first
%   step. nLow thresholds are added between max(C(1),1) and that first step,
%   evenly on a log scale, when the step is more than twice the lower value.
%
%   DUPLICATES OF THRESHOLD 0. A threshold that a one-case seed already meets
%   (0 < C_thr <= min(I0,1)) starts every sub-epidemic at the first time
%   point, exactly as threshold 0 does. Such values are removed, so the same
%   model is not ranked twice.
%
%   The result is a row, increasing, starting with 0.
%
%   See also fittingModifiedLogisticFunctionPatchABC, spatialWaveFitSettings

if nargin<4 || isempty(nLow)
    nLow = 0;
end
data = double(data(:));
n = numel(data);
if smoothfactor1>1
    smoothed = smoothdata(data,'movmean',smoothfactor1);
else
    smoothed = data;
end
C = cumsum(smoothed);
grid = linspace(C(1),C(end),n);
grid = [0 grid(1:end-1)];

if nLow>0 && numel(grid)>=3
    lower = max(grid(2),1);
    upper = grid(3);
    if upper>2*lower
        extra = exp(linspace(log(lower),log(upper),nLow+2));
        grid = [grid extra(2:end-1)];
    end
end

grid(grid>0 & grid<=min(I0,1)) = [];
onset_thrs = unique(grid(isfinite(grid) & grid>=0));
onset_thrs = onset_thrs(:)';
end
