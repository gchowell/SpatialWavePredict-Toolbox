function weights1 = getSpatialWaveAICcWeights(RMSES, topmodels1, weight_type1)
%GETSPATIALWAVEAICcWEIGHTS Compute AICc-based spatial-wave ensemble weights.
%   RMSES must have the spatial-wave layout:
%       [npatches, onset_thr, typedecline1, AICc]
%   TOPMODELS1 contains distinct, one-based row/rank indices in RMSES.
%   Output is a column vector in the same order as TOPMODELS1.
%
%   WEIGHT_TYPE1 = 0: legacy normalized inverse-AICc weights. All selected
%                    AICc values must be strictly positive.
%   WEIGHT_TYPE1 = 1: normalized Akaike relative-likelihood weights. Finite
%                    negative and zero AICc values are valid in this mode.
%
%   Only selected rows must have finite AICc. Invalid selections are rejected,
%   not silently dropped or reordered. This function does not rank models,
%   fit models, draw random samples, or change any global state.

if ~(isnumeric(RMSES) && isreal(RMSES) && ismatrix(RMSES) && ...
        ~isempty(RMSES) && size(RMSES,2) == 4)
    error('SpatialWave:InvalidRankingMatrix', ...
        ['RMSES must be a nonempty real numeric matrix with four columns: ' ...
         '[npatches onset_thr typedecline1 AICc].']);
end

if ~(isnumeric(topmodels1) && isreal(topmodels1) && ...
        isvector(topmodels1) && ~isempty(topmodels1))
    error('SpatialWave:InvalidModelSelection', ...
        'topmodels1 must be a nonempty vector of distinct valid model ranks.');
end

selectedRanks = full(double(topmodels1(:)));
if any(~isfinite(selectedRanks)) || ...
        any(selectedRanks < 1 | selectedRanks > size(RMSES,1)) || ...
        any(selectedRanks ~= fix(selectedRanks)) || ...
        numel(unique(selectedRanks)) ~= numel(selectedRanks)
    error('SpatialWave:InvalidModelSelection', ...
        'topmodels1 must contain distinct integer ranks between 1 and size(RMSES,1).');
end

if ~(isnumeric(weight_type1) && isreal(weight_type1) && ...
        isscalar(weight_type1) && any(weight_type1 == [0 1]))
    error('SpatialWave:InvalidAICcWeightType', ...
        'weight_type1 must be 0 (inverse AICc) or 1 (Akaike weights).');
end

% The fourth column is AICc; the third column is only a decline-family label.
aicc = full(double(RMSES(selectedRanks,4)));
if any(~isfinite(aicc))
    error('SpatialWave:NonfiniteAICc', ...
        'Every selected model must have a finite AICc value in RMSES(:,4).');
end

switch weight_type1
    case 0
        if any(aicc <= 0)
            error('SpatialWave:NonpositiveInverseAICc', ...
                ['Inverse-AICc weighting requires strictly positive AICc. ' ...
                 'Use weight_type1=1 for valid finite zero or negative AICc.']);
        end
        % Equivalent to 1./aicc after normalization, with no reciprocal
        % overflow for very small positive AICc values.
        relativeWeights = min(aicc) ./ aicc;

    case 1
        % Normalize relative to the best SELECTED model. At least one term
        % is exactly one, even when the selected set excludes rank 1.
        relativeWeights = exp(-0.5 * (aicc - min(aicc)));
end

weights1 = relativeWeights ./ sum(relativeWeights);
end
