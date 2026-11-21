function pii = spatialWaveIdentifiability(Phatss, LB, UB, P0, settings, B2, boundHitTol)
%SPATIALWAVEIDENTIFIABILITY Identifiability report of one ranked model.
%   pii = spatialWaveIdentifiability(Phatss, LB, UB, P0)
%   pii = spatialWaveIdentifiability(Phatss, LB, UB, P0, settings, B2, boundHitTol)
%
%   The measures the SubEpiPredict toolbox reports for each parameter
%   (subepidemicIdentifiability.m), computed from the bootstrap refits of
%   one ranked model of the spatial wave toolbox, and the correlations
%   between the parameters over the refits.
%
%   Inputs
%     Phatss    M-by-7 bootstrap refits, columns [r p a K q alpha d]
%     LB, UB    search bounds in the same layout (spatialWaveSearchSpace.m;
%               saved with a fit as bootstrapDiagnostics.searchLB/searchUB)
%     P0        the ranked fit to the observed data, 7 entries ([] = NaN)
%     settings  spatialWaveFitSettings.m ([] = the settings in force). Used:
%               PIIscale, notIdentifiedRatio, weaklyIdentifiedPII,
%               tradeoffCorrelation.
%     B2           resamples of the Monte Carlo error (default 500)
%     boundHitTol  relative tolerance of "at a bound" (default 1e-3)
%
%   Output: a struct. One entry per reported parameter (k-by-1 fields); r, p,
%   a, K and q are always listed, alpha and d only when the estimation
%   method estimates them (LB < UB).
%     name           'r', 'p', 'a', 'K', 'q', 'alpha', 'd'
%     column         column of Phatss
%     estimate       value in the ranked fit
%     median         median of the refits
%     CI_LB, CI_UB   2.5% and 97.5% quantiles of the refits
%     PII            log10(CI_UB/CI_LB), interval width in decades; NaN when CI_LB <= 0
%     PII_rangeFrac  (CI_UB - CI_LB)/(UB - LB), in [0,1]; NaN for a fixed parameter
%     PII_scale      headline measure: 1 = PII (LB > 0), 0 = PII_rangeFrac (LB <= 0)
%     PII_box        log10(UB/LB), the largest PII the bounds allow; NaN when LB <= 0
%     PII_ratio      ~1 = the interval is the search range: PII/PII_box on the
%                    log scale, PII_rangeFrac otherwise; NaN for a fixed parameter
%     fracLB, fracUB share of the refits within boundHitTol (relative) of LB, UB
%     MCSE_CI_LB, MCSE_CI_UB, MCSE_PII
%                    Monte Carlo SE of the interval limits and of the PII
%     bound_LB, bound_UB  the search bounds
%     fixed          true where LB == UB (the parameter is not estimated)
%     verdict        'fixed', 'identified', 'weakly identified', 'at a bound'
%                    or 'not identified'
%   Scalars and tables:
%     M, B2, boundHitTol, notIdentifiedRatio, weaklyIdentifiedPII,
%     tradeoffCorrelation
%     nFree          number of estimated parameters
%     nIdentified    how many of them have the verdict 'identified'
%     correlation    k-by-k rank (Spearman) correlations between the listed
%                    parameters over the refits; NaN for a fixed parameter
%                    and for one whose refits are all equal, and all NaN
%                    with fewer than 10 refits (a rank correlation of a
%                    handful of refits is +1 or -1 by chance)
%     tradeoffPairs  rows [i j rho]: pairs with |rho| >= tradeoffCorrelation
%     tradeoffs      the same as text, e.g. 'K-q (-0.97)'
%
%   THE VERDICT of an estimated parameter, in this order:
%     not identified  PII_ratio >= notIdentifiedRatio (default 0.5): the 95%
%                     interval covers at least that share of what the search
%                     bounds allow, so it reflects the bounds, not the data.
%                     (Where the PII is undefined on the log scale because
%                     CI_LB <= 0, the range fraction is used.)
%     at a bound      otherwise, more than 5% of the refits at one bound: the
%                     bound cuts the distribution of the estimate.
%     weakly identified  otherwise, an interval that is wide although the
%                     bounds do not set it: on the log scale a PII of
%                     weaklyIdentifiedPII or more (default 1: the interval
%                     spans a factor of 10 or more); on either scale a ratio
%                     of at least half of notIdentifiedRatio. With search
%                     bounds as wide as those of K (20 to 1e8), an interval
%                     from 4 000 to 400 000 covers a third of the range and
%                     would otherwise pass as identified.
%     identified      otherwise.
%
%   TRADE-OFFS. The verdict looks at one parameter at a time. Two parameters
%   that can make up for each other (K and q: a larger first sub-epidemic
%   with a faster decline; r and p) each get a wide interval although a
%   combination of them is well determined. A rank correlation near +1 or -1
%   over the refits shows that; the pairs are listed.
%
%   Reading the numbers: PII near 0 is a tight interval, PII above 1 an
%   interval wider than a factor of 10. A difference between two PII values
%   smaller than about twice MCSE_PII is bootstrap noise.
%
%   Revision SW-6 (after subepidemicIdentifiability.m of the SubEpiPredict toolbox).
%
%   See also spatialWavePII, spatialWaveBoundHitFractions,
%   spatialWaveQuantileMCSE, spatialWavePrintPII, spatialWaveWritePII

nz = 7;
Phatss = double(full(Phatss));
LB = double(full(LB(:)'));
UB = double(full(UB(:)'));
if size(Phatss,2)~=nz || numel(LB)~=nz || numel(UB)~=nz
    error('SpatialWave:ParameterDimensionMismatch', ...
        'Phatss, LB and UB must have 7 columns: [r p a K q alpha d].');
end
if nargin<4 || isempty(P0)
    P0 = nan(1,nz);
end
P0 = double(full(P0(:)'));
if numel(P0)~=nz
    error('SpatialWave:ParameterDimensionMismatch','P0 must have 7 entries.');
end
if nargin<5 || isempty(settings)
    settings = spatialWaveFitSettings();
end
scaleUser = localField(settings,'PIIscale',[]);
if isempty(scaleUser)
    scaleUser = nan(1,nz);
end
scaleUser = double(scaleUser(:)');
if numel(scaleUser)~=nz
    error('SpatialWave:PIIScale', ...
        'PIIscale must have 7 entries, one per parameter: [r p a K q alpha d].');
end
notIdentifiedRatio = localField(settings,'notIdentifiedRatio',0.5);
tradeoffCorrelation = localField(settings,'tradeoffCorrelation',0.9);
weaklyIdentifiedPII = localField(settings,'weaklyIdentifiedPII',1);
if nargin<6 || isempty(B2)
    B2 = 500;
end
if nargin<7 || isempty(boundHitTol)
    boundHitTol = 1e-3;
end

lo = quantile(Phatss,0.025,1);
hi = quantile(Phatss,0.975,1);
med = median(Phatss,1);
[PII,rangeFrac,scale,PIIbox,ratio] = spatialWavePII(lo,hi,LB,UB,scaleUser);
[fracLB,fracUB] = spatialWaveBoundHitFractions(Phatss,LB,UB,boundHitTol);
[mcseLB,mcseUB,mcsePII] = spatialWaveQuantileMCSE(Phatss,B2);

% Parameters reported: r, p, a, K, q; alpha and d when estimated.
keep = [1:5, 5+find(UB(6:7)>LB(6:7))];
names = {'r';'p';'a';'K';'q';'alpha';'d'};
pii = struct();
pii.name = names(keep);
pii.column = keep(:);
pii.estimate = P0(keep)';
pii.median = med(keep)';
pii.CI_LB = lo(keep)';
pii.CI_UB = hi(keep)';
pii.PII = PII(keep)';
pii.PII_rangeFrac = rangeFrac(keep)';
pii.PII_scale = scale(keep)';
pii.PII_box = PIIbox(keep)';
pii.PII_ratio = ratio(keep)';
pii.fracLB = fracLB(keep)';
pii.fracUB = fracUB(keep)';
pii.MCSE_CI_LB = mcseLB(keep)';
pii.MCSE_CI_UB = mcseUB(keep)';
pii.MCSE_PII = mcsePII(keep)';
pii.bound_LB = LB(keep)';
pii.bound_UB = UB(keep)';
pii.fixed = (LB(keep)==UB(keep))';
pii.M = size(Phatss,1);
pii.B2 = B2;
pii.boundHitTol = boundHitTol;
pii.notIdentifiedRatio = notIdentifiedRatio;
pii.tradeoffCorrelation = tradeoffCorrelation;
pii.weaklyIdentifiedPII = weaklyIdentifiedPII;

% ---- verdict
nk = numel(keep);
pii.verdict = cell(nk,1);
for k = 1:nk
    measure = pii.PII_ratio(k);
    if isnan(measure)
        measure = pii.PII_rangeFrac(k);   % log scale asked for, CI_LB <= 0
    end
    if pii.fixed(k)
        pii.verdict{k} = 'fixed';
    elseif measure>=notIdentifiedRatio
        pii.verdict{k} = 'not identified';
    elseif pii.fracLB(k)>0.05 || pii.fracUB(k)>0.05
        pii.verdict{k} = 'at a bound';
    elseif (pii.PII_scale(k)==1 && pii.PII(k)>=weaklyIdentifiedPII) || measure>=notIdentifiedRatio/2
        pii.verdict{k} = 'weakly identified';
    else
        pii.verdict{k} = 'identified';
    end
end
pii.nFree = sum(~pii.fixed);
pii.nIdentified = sum(strcmp(pii.verdict,'identified'));

% ---- rank correlations between the parameters over the refits
ranks = nan(size(Phatss,1),nk);
for k = 1:nk
    x = Phatss(:,keep(k));
    if ~pii.fixed(k) && all(isfinite(x)) && max(x)>min(x)
        ranks(:,k) = localRanks(x);
    end
end
pii.correlation = nan(nk,nk);
pii.tradeoffPairs = zeros(0,3);
pii.tradeoffs = {};
for i = 1:nk
    for j = i:nk
        if all(isfinite(ranks(:,i))) && all(isfinite(ranks(:,j))) && size(ranks,1)>=10
            a = ranks(:,i)-mean(ranks(:,i));
            b = ranks(:,j)-mean(ranks(:,j));
            rho = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
            pii.correlation(i,j) = rho;
            pii.correlation(j,i) = rho;
            if j>i && abs(rho)>=tradeoffCorrelation
                pii.tradeoffPairs(end+1,:) = [i j rho]; %#ok<AGROW>
                pii.tradeoffs{end+1,1} = sprintf('%s-%s (%+.2f)',pii.name{i},pii.name{j},rho); %#ok<AGROW>
            end
        end
    end
end
end

% -----------------------------------------------------------------------------
function v = localField(S,name,default)
if isstruct(S) && isfield(S,name)
    v = S.(name);
else
    v = default;
end
end

% -----------------------------------------------------------------------------
function r = localRanks(x)
% Ranks 1..M with ties given their average rank.
[sorted,order] = sort(x(:));
M = numel(sorted);
r = zeros(M,1);
first = 1;
while first<=M
    last = first;
    while last<M && sorted(last+1)==sorted(first)
        last = last+1;
    end
    r(order(first:last)) = (first+last)/2;
    first = last+1;
end
end
