function assertSpatialWaveFitIdentity(s,rank1,ranking)
%ASSERTSPATIALWAVEFITIDENTITY Reject old/mixed model or forecast rank files.
% S may be a loaded MAT struct or just a struct with fitRunId and fitResult.
if ~isstruct(s) || ~all(isfield(s,{'fitRunId','fitResult'})) || ...
        ~strcmp(s.fitRunId,ranking.fitRunId)
    error('SpatialWave:StaleFinalFitResults', ...
        ['This model/forecast file is legacy or belongs to a different fit run. ' ...
         'Rerun fitting and regenerate all dependent forecasts before ensembling.']);
end
if ~isscalar(rank1) || rank1<1 || rank1>ranking.numModels || rank1~=fix(rank1)
    error('SpatialWave:InvalidModelSelection','Requested rank is outside the accepted-fit ranking.');
end
f=s.fitResult; reference=ranking.fits{rank1};
fields={'rank','parameters','modelSpec','AICc','objectiveValue','bestfit','target','data'};
for j=1:numel(fields)
    key=fields{j};
    if ~isfield(f,key) || ~isequaln(f.(key),reference.(key))
        error('SpatialWave:StaleFinalFitResults', ...
            'The saved model/forecast does not match accepted final rank %d (%s).',rank1,key);
    end
end
end
