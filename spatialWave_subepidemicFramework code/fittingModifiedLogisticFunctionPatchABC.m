% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [RMSES,PS,npatches,onset_thr,typedecline1,P]=fittingModifiedLogisticFunctionPatchABC( ...
    datafilename1,data1,DT,epidemic_period,M,flagX,typedecline2,numstartpoints)
%FITTINGMODIFIEDLOGISTICFUNCTIONPATCHABC Search over the candidate spatial wave models.
%   Fits every candidate (number of sub-epidemics, onset threshold, decline
%   family) and ranks the fits by AICc.
%
%   RMSES  one row per fitted candidate, best first:
%          [active sub-epidemics, C_thr, decline family, AICc]
%   PS     the fitted parameters of each row, [r p a K q alpha d]
%
%   Revision SW-2. Every candidate is fitted by fitSpatialWaveModel.m: only
%   the free parameters are optimized, each mapped to [0,1] (r, K and alpha on
%   a log scale); the objective gets its model from an explicit context, not
%   from global variables; every local minimum is evaluated again before it
%   is accepted; and each fit returns diagnostics. What else changed here:
%     * the fit of each threshold also starts from the solution of the
%       previous threshold (setting continuationStart);
%     * thresholds are added at the low end of the grid
%       (spatialWaveThresholdGrid.m, setting thresholdLowPoints);
%     * the single sub-epidemic candidate is fitted once, not once per decline
%       family: the family has no effect on one sub-epidemic, and the two
%       fits were the same model ranked twice;
%     * a candidate without a converged, valid fit is left out and reported;
%     * count likelihoods are fitted to the raw counts (spatialWaveFitSeries.m).
%   Settings: spatialWaveFitSettings.m and options_fit_settings_sw.m.
%
%   The global variables below are still written, for scripts that read
%   them; the fit itself no longer reads any of them.

global flag1 method1 timevect ydata yfit
global I0 npatches onset_thr typedecline1

flag1=flagX;
close all

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed
global onset_fixed

global dist1
global factor1

global smoothfactor1
global calibrationperiod1

global LBe UBe

fitSettings=spatialWaveFitSettings();
odeOptions=spatialWaveOdeOptions();

% <=============================================================================================>
% <========= Bounds of alpha and d (kept in the globals for older scripts) =====================>
% <=============================================================================================>
switch method1
    case {0,1,2,6}
        LBe=[0 0];
        UBe=[0 0];
    case {3,4}
        LBe=[10^-8 0];
        UBe=[10^5 0];
    case 5
        LBe=[10^-8 0.2];
        UBe=[10^5 10^2];
end

% <==============================================================================>
% <============ Load data and proceed to parameter estimation ===================>
% <==============================================================================>
data1=data1(epidemic_period,:);
if data1(1,2)==0
    data1=data1(2:end,:);
end
data = data1(:,2);

% <==============================================================================>
% <============ Set time vector (timevect) and initial condition (I0) ===========>
% <==============================================================================>
timevect = (data1(:,1));
I0 = data(1); % initial condition

% <==============================================================================>
% <============ Series the objective is fitted to ===============================>
% <==============================================================================>
% Least squares and absolute deviations: the moving average of span
% smoothfactor1, as before. Count likelihoods (methods 1, 3, 4, 5): the raw
% counts; smoothfactor1 then only shapes the grid of onset thresholds.
[ydata,fitSeriesInfo]=spatialWaveFitSeries(data,smoothfactor1,method1,fitSettings);

% Revision SW-5: day-of-week factors (spatialWaveWeekdayFactors.m).
% weekdayCurve is [] when no adjustment applies, and the search is then as
% before. With the adjustment the model describes the series without the
% weekly cycle, so its initial condition is the first observation divided by
% its factor.
weekday=spatialWaveWeekdayFactors(data,fitSettings,fitSeriesInfo);
weekdayCurve=spatialWaveWeekdayCurve(weekday,length(data));
if ~isempty(weekdayCurve)
    I0=I0/weekdayCurve(1);
    fprintf('Day-of-week factors, from the day of the first observation: %s\n',mat2str(weekday.byPosition',3));
end
if fitSeriesInfo.smoothedCountLikelihood
    warning('SpatialWave:SmoothedCountLikelihood', ...
        ['method1 = %d fits a count likelihood to a %g-point moving average. The fitted dispersion ', ...
         'will be too small and prediction intervals too narrow for the raw counts.'],method1,smoothfactor1);
elseif any(method1==[1 3 4 5]) && smoothfactor1>1
    fprintf(['method1 = %d is fitted to the raw counts. smoothfactor1 = %g is used only for the grid of ', ...
        'onset thresholds.\n'],method1,smoothfactor1);
end

% <==============================================================================>
% <============ Search space, seed and start points =============================>
% <==============================================================================>
% The parameters are shared by all sub-epidemics, so the bounds, the seed and
% the start points are the same for every candidate. Rows of starts_base are
% in original units; the Latin hypercube is drawn in unit coordinates, so r,
% K and alpha are sampled log-uniformly.
[LB,UB,isLog,z0]=spatialWaveSearchSpace(data,flag1,method1,fitSettings);
P0=z0; % earlier name of the seed, kept in the saved workspace
unitMap=spatialWaveUnitMap('make',struct('z0',z0,'LB',LB,'UB',UB,'isLog',isLog));
nStarts=max(0,floor(numstartpoints));
starts_base=spatialWaveStartPoints(unitMap,nStarts,0,z0);

% <==============================================================================>
% <================= Onset thresholds and numbers of sub-epidemics ==============>
% <==============================================================================>
onset_thrs2=spatialWaveThresholdGrid(data,smoothfactor1,I0,fitSettings.thresholdLowPoints);

if onset_fixed==1
    npatchess = 1:1:npatches_fixed;
else
    npatchess = unique([1 npatches_fixed]);
end

RMSES = zeros(0,4);
PS    = zeros(0,7);

% One row per fitted candidate, in search order, rejected ones included.
% candidateFitInfo{k} is the full fitInfo structure of row k.
fitDiagnosticsColumns={'configured','onset_thr','typedecline','success','exitflag','nActive', ...
    'startAgreement','gradientInfNorm','polishRun','polishImproved','nLocalRuns','nConvergedRuns','nAtBound','fval'};
fitDiagnostics=zeros(0,numel(fitDiagnosticsColumns));
candidateFitInfo={};
rejectedCandidates=0;
candidatesTooLarge=zeros(0,3); % [configured sub-epidemics, parameters, data points]

% <================================================================================================>
% <==== Evaluate AICc across models with different number of subepidemics and C_thr values ========>
% <================================================================================================>
for npatches2 = npatchess
    % Keep the configured search dimension separate from activated counts.
    configuredNpatches = npatches2;

    % AICc needs more data points than parameters plus one (getAICc.m). A
    % candidate that is too large for this series is not fitted.
    kCandidate=get_nparams(method1,dist1,configuredNpatches,flag1,1,onset_fixed);
    if length(ydata)<=kCandidate+1
        candidatesTooLarge(end+1,:)=[configuredNpatches kCandidate length(ydata)]; %#ok<AGROW>
        warning('SpatialWave:CandidateTooLarge', ...
            ['Models with %d sub-epidemic(s) have %d parameters, too many for the %d data points ', ...
             '(AICc needs more than %d). They are left out of the search.'], ...
            configuredNpatches,kCandidate,length(ydata),kCandidate+1);
        continue
    end

    if (onset_fixed==1 || configuredNpatches==1)
        onset_thrs = 0;
    else
        onset_thrs = onset_thrs2;
    end

    % One sub-epidemic has one size: the decline family has no effect on it.
    if configuredNpatches==1
        families = typedecline2(1);
    else
        families = typedecline2;
    end

    for j = 1:length(families)

        % Accepted solution of the previous threshold of this family
        % (continuation start); empty until one fit has been accepted.
        previousP=[];

        for onset_thr = onset_thrs
            % Each threshold/family candidate is fitted with its own configured
            % dimension, irrespective of how many components activated before.
            npatches = configuredNpatches;
            typedecline1 = families(j);

            % Explicit fit context: the objective reads nothing from globals.
            context=struct('timevect',timevect,'ydata',ydata,'I0',I0, ...
                'npatches',configuredNpatches,'onset_fixed',onset_fixed,'onset_thr',onset_thr, ...
                'flag1',flag1,'typedecline1',typedecline1,'method1',method1,'odeOptions',odeOptions, ...
                'weekdayFactor',weekdayCurve);

            % The solution of the previous threshold is an excellent start for
            % this one. It is one more start point, never a replacement. A fit
            % without one gets firstThresholdStartFactor times the start points.
            starts=starts_base;
            nStartsThisFit=nStarts;
            if fitSettings.continuationStart && ~isempty(previousP)
                starts=[starts_base;previousP];
            elseif isempty(previousP) && fitSettings.firstThresholdStartFactor>1 && nStarts>0
                nExtra=ceil((fitSettings.firstThresholdStartFactor-1)*nStarts);
                extraStarts=spatialWaveStartPoints(unitMap,nExtra,0,z0);
                starts=[extraStarts(1:end-1,:);starts_base]; % the seed is already the last row of starts_base
                nStartsThisFit=nStarts+nExtra;
            end

            [P,fval,fitInfo]=fitSpatialWaveModel(context,z0,LB,UB,isLog,starts,nStartsThisFit,fitSettings);

            candidateFitInfo{end+1,1}=fitInfo; %#ok<AGROW>
            fitDiagnostics(end+1,:)=[configuredNpatches onset_thr typedecline1 double(fitInfo.success) ...
                fitInfo.exitflag fitInfo.nActive fitInfo.startAgreement fitInfo.gradientInfNorm ...
                double(fitInfo.polish.run) double(fitInfo.polish.improved) fitInfo.nLocalRuns ...
                fitInfo.nConvergedRuns numel(fitInfo.atLowerBound)+numel(fitInfo.atUpperBound) fval]; %#ok<AGROW>

            % A candidate with no converged, valid local minimum is recorded
            % and left out. It is never replaced by its seed.
            if ~fitInfo.success
                rejectedCandidates=rejectedCandidates+1;
                continue
            end
            previousP=P;

            % Preserve the existing activated-count reporting/AICc policy,
            % without changing the configured dimension used by later fits.
            activatedNpatches = min(configuredNpatches,fitInfo.nActive);

            AICc = getAICc(method1,dist1,activatedNpatches,flag1(1),1,fval,length(ydata),onset_fixed);

            RMSES(end+1,:) = [activatedNpatches onset_thr typedecline1 AICc]; %#ok<AGROW>
            PS(end+1,:)    = P; %#ok<AGROW>
        end % onset_thr
    end % typedecline1
end % npatches2

if rejectedCandidates>0
    warning('SpatialWave:CandidatesRejected', ...
        ['%d of %d candidate model(s) had no converged, valid fit and were left out. ', ...
         'fitDiagnostics in the ABC file lists them.'],rejectedCandidates,size(fitDiagnostics,1));
end
if isempty(RMSES) && ~isempty(candidatesTooLarge) && isempty(fitDiagnostics)
    error('SpatialWave:TooFewDataPoints', ...
        ['No model can be ranked: the series has %d data points and the smallest model has %d parameters. ', ...
         'AICc needs more data points than parameters plus one.'],length(ydata),min(candidatesTooLarge(:,2)));
end
if isempty(RMSES)
    error('SpatialWave:NoValidCandidates','No candidate model has a converged, valid fit.');
end

% <=============================================================================================>
% <======================== Sort the results by AICc (lowest to highest) =======================>
% <=============================================================================================>
[RMSES,index1] = sortrows(RMSES,[4 1]);
PS = PS(index1,:);

[~, index1] = min(RMSES(:,4));

npatches     = RMSES(index1,1);
onset_thr    = RMSES(index1,2);
typedecline1 = RMSES(index1,3);
AICmin       = RMSES(index1,4);
relativelik_i = exp((AICmin - RMSES(:,4))/2);

% <=============================================================================================>
% <============================ Get the best fit resultls ======================================>
% <=============================================================================================>
P = PS(index1,:);

r_hat=P(1,1);
p_hat=P(1,2);
a_hat=P(1,3);
K_hat=P(1,4);
q_hat=P(1,5);

alpha_hat=P(1,end-1);
d_hat=P(1,end);

if method1==3
    dist1=3;  factor1=alpha_hat; % VAR=mean+alpha*mean;
elseif method1==4
    dist1=4;  factor1=alpha_hat; % VAR=mean+alpha*mean^2;
elseif method1==5
    dist1=5;  factor1=alpha_hat; % VAR=mean+alpha*mean^2;
end

% The search fitted this row with its configured number of sub-epidemics;
% RMSES holds the number that became active (see spatialWavePatchCap.m).
patchCap=spatialWavePatchCap(npatches,onset_thr,onset_fixed,npatches_fixed);
[~,x,totinc,patchIncidence]=simulateSpatialWave(timevect,I0,patchCap,onset_fixed, ...
             onset_thr,flag1,typedecline1,r_hat,p_hat,a_hat,K_hat,q_hat);
if ~isempty(weekdayCurve) % revision SW-5: the curve on the scale of the data
    totinc=totinc(:).*weekdayCurve;
    patchIncidence=bsxfun(@times,patchIncidence,weekdayCurve);
end

% <=============================================================================================>
% <================================= Plot best model fit =======================================>
% <=============================================================================================>
figure(100)
for j=1:npatches
    plot(timevect,patchIncidence(:,j)); hold on
end
bestfit=totinc;

plot(timevect,totinc,'r'); hold on
plot(timevect,data,'bo');
xlabel('Time (days)'); ylabel('Incidence');
title('best fit'); [npatches onset_thr];

% overdispersion estimate if needed
if (method1==0 && dist1==2)
    binsize1=7;
    [ratios,~]=getMeanVarianceRatio(data,binsize1,2);
    index1=find(ratios(:,1)>0);
    factor1=mean(ratios(index1,1));
end

'ABC estimates:'
'[npatches onset_thr typedecline1]'
[npatches onset_thr typedecline1]

% <=============================================================================================>
% <===================================  Save the results  ======================================>
% <=============================================================================================>
save(strcat('./output/ABC-original-npatchesfixed-',num2str(npatches_fixed), ...
    '-onsetfixed-',num2str(onset_fixed),'-typedecline-',num2str(sum(typedecline2)), ...
    '-smoothing-',num2str(smoothfactor1),'-',datafilename1(1:end-4), ...
    '-flag1-',num2str(flag1(1)),'-method-',num2str(method1), ...
    '-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'.mat'),'-mat')

end
