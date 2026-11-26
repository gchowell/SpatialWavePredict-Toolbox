function spatialWavePrintPII(pii, heading)
%SPATIALWAVEPRINTPII Print the identifiability report of one ranked model.
%   spatialWavePrintPII(pii, heading) prints the struct returned by
%   spatialWaveIdentifiability.m.
%
%   Columns: estimate and 95% interval; PII = log10(UB/LB) of the interval;
%   range = share of the search range the interval covers; scale = headline
%   measure (log: PII, range: range fraction); box = largest PII the bounds
%   allow; ratio = PII/box (log) or range fraction, ~1 when the interval is
%   the search range; atLB/atUB = share of refits at a bound; MCSE = Monte
%   Carlo SE of the PII; then the verdict. Below the table: the pairs of
%   parameters that trade off against each other.
%
%   Revision SW-6 (after subepidemicPrintPII.m of the SubEpiPredict toolbox).
%
%   See also spatialWaveIdentifiability, spatialWaveWritePII

if nargin<2 || isempty(heading)
    heading = 'Parameter identifiability (PII)';
end
fprintf('\n%s (%d bootstrap refits):\n',heading,pii.M);
fprintf('  %-7s %10s %10s %10s %7s %7s %6s %7s %7s %6s %6s %7s  %s\n', ...
    'param','estimate','CI LB','CI UB','PII','range','scale','box','ratio','atLB','atUB','MCSE','verdict');
for k = 1:numel(pii.name)
    if pii.fixed(k)
        fprintf('  %-7s %10.4g %s\n',pii.name{k},pii.estimate(k),'  fixed');
        continue
    end
    if pii.PII_scale(k)==1
        scaleName = 'log';
    else
        scaleName = 'range';
    end
    fprintf('  %-7s %10.4g %10.4g %10.4g %7.3f %7.3f %6s %7.3f %7.3f %6.3f %6.3f %7.3f  %s\n', ...
        pii.name{k},pii.estimate(k),pii.CI_LB(k),pii.CI_UB(k),pii.PII(k),pii.PII_rangeFrac(k), ...
        scaleName,pii.PII_box(k),pii.PII_ratio(k),pii.fracLB(k),pii.fracUB(k),pii.MCSE_PII(k),pii.verdict{k});
end
fprintf('  PII = log10(CI UB/CI LB); NaN when CI LB <= 0 (read the range fraction). ratio near 1 or atLB/atUB above 0.05: the bounds set the interval.\n');
fprintf('  Identified: %d of %d estimated parameter(s). Not identified: ratio of %.2g or more. Weakly: PII of %.2g or more, or ratio of %.2g or more.\n', ...
    pii.nIdentified,pii.nFree,pii.notIdentifiedRatio,pii.weaklyIdentifiedPII,pii.notIdentifiedRatio/2);
if pii.M<50
    fprintf('  Only %d refits: the interval limits are rough (see MCSE); B of 100 or more is advisable for this report.\n',pii.M);
end
if pii.M<10
    fprintf('  Rank correlations between the parameters need at least 10 refits.\n\n');
elseif isempty(pii.tradeoffs)
    fprintf('  No pair of parameters with a rank correlation of %.2g or more (in size) over the refits.\n\n',pii.tradeoffCorrelation);
else
    fprintf('  Parameters that trade off against each other (rank correlation over the refits): %s\n\n',strjoin(pii.tradeoffs',', '));
end
end
