function pii = spatialWaveSavedIdentifiability(Phatss, P0, bootstrapDiagnostics)
%SPATIALWAVESAVEDIDENTIFIABILITY Identifiability report of a saved fit.
%   pii = spatialWaveSavedIdentifiability(Phatss, P0, bootstrapDiagnostics)
%   computes spatialWaveIdentifiability.m from what the file of a ranked
%   model holds: the bootstrap refits, the ranked fit and the search bounds
%   (bootstrapDiagnostics.searchLB/searchUB, saved since revision SW-2).
%   The report settings in force now are used (PIIscale,
%   notIdentifiedRatio, weaklyIdentifiedPII, tradeoffCorrelation), so no fit has to be repeated
%   to change them.
%
%   Returns [] for a fit that does not hold its search bounds (made before
%   revision SW-2) or whose arrays do not have 7 columns.
%
%   Revision SW-6.
pii = [];
if ~(isstruct(bootstrapDiagnostics) && isfield(bootstrapDiagnostics,'searchLB') && isfield(bootstrapDiagnostics,'searchUB'))
    return
end
if ~(isnumeric(Phatss) && size(Phatss,2)==7 && size(Phatss,1)>=1 && numel(P0)==7 && ...
        numel(bootstrapDiagnostics.searchLB)==7 && numel(bootstrapDiagnostics.searchUB)==7)
    return
end
pii = spatialWaveIdentifiability(Phatss,bootstrapDiagnostics.searchLB,bootstrapDiagnostics.searchUB,P0);
end
