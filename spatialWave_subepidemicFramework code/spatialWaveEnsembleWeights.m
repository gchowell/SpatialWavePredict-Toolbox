function weights = spatialWaveEnsembleWeights(weight_type1,RMSES,topmodels1,WISC,WISF)
%SPATIALWAVEENSEMBLEWEIGHTS Weights of the member models of an ensemble.
%   weights = spatialWaveEnsembleWeights(weight_type1,RMSES,topmodels1,WISC,WISF)
%
%   RMSES is the ranking matrix of the spatial-wave fits,
%       [npatches onset_thr typedecline1 AICc],
%   and topmodels1 holds the ranks (rows of RMSES) of the member models.
%   WISC and WISF hold the weighted interval scores of the members
%   (calibration; earlier forecasts), one per member and in the order of
%   topmodels1. They are read only for codes 2 and 3; pass [] otherwise.
%   The result is a column of nonnegative weights that sum to 1, one per
%   member, in the order of topmodels1.
%
%   weight_type1
%     -1  equal weights
%      0  Akaike weights (see below)
%      1  Akaike weights: exp(-(AICc - min(AICc))/2), normalized
%      2  proportional to 1/WISC
%      3  proportional to 1/WISF
%
%   Code 0 used to be (1/AICc)/sum(1/AICc). That is not a weight. AICc is
%   defined up to an additive constant, so only AICc differences carry
%   information: models whose AICc are large and close get almost equal
%   weights whatever their differences, and the code stopped with an error
%   whenever an AICc was zero or negative, which a likelihood-based fit can
%   give. Code 0 now returns Akaike weights and warns
%   (SpatialWave:WeightCodeZero); use code 1 for the same weights without
%   the warning, or -1 for equal weights. The inverse-AICc numbers are still
%   available from getSpatialWaveAICcWeights(RMSES,topmodels1,0).
%
%   Scores of codes 2 and 3 must be finite and positive; anything else is an
%   error, because 1/score would be infinite or would change sign. They are
%   zero when scoring is off (getperformance = 0).
%
%   Added in revision SW-10.
%
%   See also getensemblesubepidemics, getSpatialWaveAICcWeights, options_forecast

if nargin<4
    WISC=[];
end
if nargin<5
    WISF=[];
end

if ~(isnumeric(weight_type1) && isreal(weight_type1) && isscalar(weight_type1) && ...
        any(weight_type1==[-1 0 1 2 3]))
    error('SpatialWave:EnsembleWeights', ...
        'weight_type1 must be -1, 0, 1, 2 or 3 (see options_forecast.m).');
end
if ~(isnumeric(topmodels1) && isreal(topmodels1) && isvector(topmodels1) && ~isempty(topmodels1))
    error('SpatialWave:InvalidModelSelection', ...
        'topmodels1 must be a nonempty vector of distinct valid model ranks.');
end
k = numel(topmodels1);

switch weight_type1
    case -1
        weights = ones(k,1)/k;
    case {0,1}
        if weight_type1==0
            warning('SpatialWave:WeightCodeZero', ...
                ['weight_type1 = 0 used 1/AICc, which is not a weight (it depends on the additive ', ...
                 'constant of AICc and fails when AICc is zero or negative). It now gives Akaike ', ...
                 'weights, as weight_type1 = 1 does.']);
        end
        % AICc is read from column 4 of RMSES; the selection is validated there.
        weights = getSpatialWaveAICcWeights(RMSES,topmodels1,1);
    case 2
        weights = localInverseScore(WISC,k,'WISC',2);
    case 3
        weights = localInverseScore(WISF,k,'WISF',3);
end

weights = weights(:);
if ~(all(isfinite(weights)) && all(weights>=0) && abs(sum(weights)-1)<1e-12)
    error('SpatialWave:EnsembleWeights', ...
        'The ensemble weights are not finite, nonnegative numbers that sum to 1.');
end
end

% -----------------------------------------------------------------------------
function weights = localInverseScore(score,k,name,code)
if ~(isnumeric(score) && isreal(score) && numel(score)==k)
    error('SpatialWave:EnsembleWeights', ...
        'weight_type1 = %d needs one %s value per member model (%d).',code,name,k);
end
score = full(double(score(:)));
if ~(all(isfinite(score)) && all(score>0))
    error('SpatialWave:EnsembleWeights', ...
        ['weight_type1 = %d needs finite, positive %s values; received %s. ', ...
         'They are zero when getperformance = 0.'],code,name,mat2str(score.',6));
end
weights = (1./score)/sum(1./score);
end
