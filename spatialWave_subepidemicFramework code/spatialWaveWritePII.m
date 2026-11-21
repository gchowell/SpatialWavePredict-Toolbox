function spatialWaveWritePII(filename, pii, ranks)
%SPATIALWAVEWRITEPII Write identifiability reports to a CSV file.
%   spatialWaveWritePII(filename, pii) writes one report
%   (spatialWaveIdentifiability.m); with a cell array of reports and RANKS,
%   one block per ranked model with its rank in the first column.
%
%   One row per parameter: the columns of the printed table at full
%   precision, the search bounds, whether the parameter is fixed, the
%   verdict, and the parameter it is most strongly correlated with over the
%   refits (rank correlation).
%
%   Revision SW-6 (after subepidemicWritePII.m of the SubEpiPredict toolbox).
%
%   See also spatialWaveIdentifiability, spatialWavePrintPII

if isstruct(pii)
    pii = {pii};
end
withRank = nargin>=3 && ~isempty(ranks);
if withRank && numel(ranks)~=numel(pii)
    error('SpatialWave:PIIRanks','ranks must have one entry per report.');
end
cols = {'estimate','median','CI_LB','CI_UB','PII','PII_rangeFrac','PII_scale','PII_box', ...
    'PII_ratio','fracLB','fracUB','MCSE_CI_LB','MCSE_CI_UB','MCSE_PII','bound_LB','bound_UB'};
fid = fopen(filename,'w');
if fid<0
    error('SpatialWave:PIIFile','The file %s could not be opened for writing.',filename);
end
if withRank
    fprintf(fid,'ranked_model,');
end
fprintf(fid,'parameter');
fprintf(fid,',%s',cols{:});
fprintf(fid,',fixed,verdict,most_correlated_with,rank_correlation\n');
for m = 1:numel(pii)
    S = pii{m};
    for k = 1:numel(S.name)
        if withRank
            fprintf(fid,'%d,',ranks(m));
        end
        fprintf(fid,'%s',S.name{k});
        for c = 1:numel(cols)
            fprintf(fid,',%s',localNumber(S.(cols{c})(k)));
        end
        other = abs(S.correlation(k,:));
        other(k) = NaN;
        if all(isnan(other))
            partner = ''; rho = NaN;
        else
            [~,j] = max(other);          % max ignores NaN
            partner = S.name{j}; rho = S.correlation(k,j);
        end
        fprintf(fid,',%d,%s,%s,%s\n',S.fixed(k),S.verdict{k},partner,localNumber(rho));
    end
end
fclose(fid);
end

% -----------------------------------------------------------------------------
function s = localNumber(x)
if isnan(x)
    s = 'NaN';
elseif isinf(x)
    if x>0, s = 'Inf'; else, s = '-Inf'; end
else
    s = sprintf('%.10g',x);
end
end
