% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [Phatss,npatches,onset_thr,typedecline1,curves,bestfit,data1,P0,AICc_best,RelLik_best,factor1,d,fitResult,bootstrapDiagnostics]=fittingModifiedLogisticFunctionPatchMultiple(RMSES,relativelik_i,PS,data1,DT,epidemic_period,M,flagX,numstartpoints,rank1)
%FITTINGMODIFIEDLOGISTICFUNCTIONPATCHMULTIPLE Accepted fit and bootstrap of one ranked model.
%   Takes row RANK1 of the search ranking (RMSES, PS), scores its curve, and
%   derives parameter uncertainty by parametric bootstrap: M synthetic
%   datasets are drawn around the fitted curve and each is refitted.
%
%   Revision SW-2.
%     * The ranked fit is bootstrapped as it stands. Before, its parameters
%       were optimized once more from a single start; when that run moved,
%       the bootstrap was centred on a fit that no longer had the AICc it was
%       ranked with. The search now returns converged fits
%       (fitSpatialWaveModel.m), so the ranked vector is used unchanged and
%       its score is checked.
%     * Every refit goes through fitSpatialWaveModel.m with an explicit fit
%       context: scaled search space, no global variables. The model
%       (number of sub-epidemics, threshold, decline family) is passed in
%       the context and cannot be changed by another fit.
%     * The synthetic datasets and the start points of every refit are drawn
%       first, on the MATLAB client; the refits can then run on a parallel
%       pool (setting useParallel) and give the same numbers as serially.
%   bootstrapDiagnostics (new, last output) records what happened to each
%   refit.
%
%   The global variables below are still written for scripts that read them.

global flag1 timevect ydata yfit

% (The fit no longer reads these; the selected model is passed in the context.)
global I0 npatches onset_thr typedecline1

flag1=flagX;

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed

global onset_fixed

global subepicurves

global method1
global dist1
global factor1

global smoothfactor1

global LBe UBe


close all

fitSettings=spatialWaveFitSettings();
odeOptions=spatialWaveOdeOptions();

% <=============================================================================================>
% <========== Load previously inferred models according to AICc (best to worst fits) ===========>
% <=============================================================================================>

% remove repeated rows
[RMSES,index1]=unique(RMSES,'rows','stable');
PS=PS(index1,:);
relativelik_i=relativelik_i(index1);

index1=rank1;

npatches=RMSES(index1,1);

onset_thr=RMSES(index1,2);

typedecline1=RMSES(index1,3);

if ~isscalar(typedecline1) || ~isfinite(typedecline1) || ~ismember(typedecline1,[1 2])
    error('SpatialWave:InvalidDeclineFamily', ...
        ['The selected decline family must be 1 (exponential) or ' ...
         '2 (power-law). Received %g.'], typedecline1);
end

% Number of sub-epidemics the search fitted this candidate with. Column 1 of
% RMSES is the number that became ACTIVE, which can be smaller. The bootstrap
% and the forecast use the fitted number, so a bootstrap dataset or a
% forecast can start a sub-epidemic that the fitted threshold allows
% (revision SW-1; see spatialWavePatchCap.m).
nActiveSearch=npatches;
patchCap=spatialWavePatchCap(nActiveSearch,onset_thr,onset_fixed,npatches_fixed);
npatches=patchCap;

searchAICc=RMSES(index1,4);
searchRelativeLikelihood=relativelik_i(index1);
% A single fit has no valid final relative likelihood until all selected
% accepted fits have been rescored. The runner fills this after reranking.
RelLik_best=NaN;

P0=full(double(PS(index1,:)));
P0_search=P0;

% <=============================================================================================>
% <================================ Data and fit context ======================================>
% <=============================================================================================>

data1=data1(epidemic_period,:);

I0=data1(1,2); % initial condition

if I0==0

    data1=data1(2:end,:);

end

data=data1(:,2);

timevect=data1(:,1);

I0=data(1); % initial condition

% Same bounds, coordinate scales and fitted series as the search.
[LB,UB,isLog]=spatialWaveSearchSpace(data,flag1,method1,fitSettings);
[ydata,fitSeriesInfo]=spatialWaveFitSeries(data,smoothfactor1,method1,fitSettings);

% Revision SW-5: day-of-week factors, as in the search. weekdayCurve is []
% when no adjustment applies, and everything below is then as before. With
% the adjustment the model describes the series without the weekly cycle, so
% its initial condition is the first observation divided by its factor.
weekday=spatialWaveWeekdayFactors(data,fitSettings,fitSeriesInfo);
weekdayCurve=spatialWaveWeekdayCurve(weekday,length(data));
if ~isempty(weekdayCurve)
    I0=I0/weekdayCurve(1);
end

% Explicit fit context: the objective reads nothing from globals. The decline
% family of the selected model, RMSES(index1,3), is part of it.
context=struct('timevect',timevect,'ydata',ydata,'I0',I0, ...
    'npatches',patchCap,'onset_fixed',onset_fixed,'onset_thr',onset_thr, ...
    'flag1',flag1,'typedecline1',typedecline1,'method1',method1,'odeOptions',odeOptions, ...
    'weekdayFactor',weekdayCurve);

% <=============================================================================================>
% <===== The ranked fit, as it stands =========================================================>
% <=============================================================================================>

P=P0;
[fval,rankedDetail]=spatialWaveObjective(P,context);
if ~rankedDetail.success
    error('SpatialWave:RankedFitInvalid', ...
        'The ranked fit %d cannot be scored: %s',rank1,rankedDetail.message);
end

r_hat=P(1);
p_hat=P(2);
a_hat=P(3);
K_hat=P(4);
q_hat=P(5);

alpha_hat=P(end-1);
d_hat=P(end);

% The accepted curve comes from the simulator the objective used, with the
% same number of sub-epidemics, so it is the curve that was optimized.
[~,x,totinc,patchIncidence,activation]=simulateSpatialWave(timevect,I0,patchCap,onset_fixed, ...
    onset_thr,flag1,typedecline1,r_hat,p_hat,a_hat,K_hat,q_hat,odeOptions);
if ~isempty(weekdayCurve) % revision SW-5: the curve on the scale of the data
    totinc=totinc(:).*weekdayCurve;
    patchIncidence=bsxfun(@times,patchIncidence,weekdayCurve);
end
nActiveFit=activation.nActive; % sub-epidemics active by the end of the calibration period

figure(10)

for j=1:nActiveFit

    plot(timevect,patchIncidence(:,j))
    hold on

end

bestfit=totinc;

% Freeze a matched accepted result BEFORE the bootstrap. Score the actual
% displayed/saved curve with the arithmetic of the objective.
% The AICc counts the sub-epidemics that are active in the accepted curve
% (the policy of the search); patchCap is stored beside it.
rankedOutput=struct('message','Ranked fit of the search, used unchanged (revision SW-2).');
fitResult=makeSpatialWaveAcceptedFit(P,P0_search,bestfit,x(:,1:nActiveFit),data1,ydata, ...
    method1,dist1,nActiveFit,flag1,onset_thr,typedecline1,onset_fixed, ...
    fval,1,rankedOutput,rank1,searchAICc,searchRelativeLikelihood,patchCap);
AICc_best=fitResult.AICc;
numparams=fitResult.numParameters;

% normal distribution of the error structure
if method1==0 & dist1==0

    var1=sum((bestfit-data).^2)./(length(bestfit)-numparams); % last revised: 16 Oct 2022
    factor1=sqrt(var1);

elseif method1==6 & dist1==6

    var1=sum((bestfit-data).^2)./(length(bestfit)-numparams); % last revised: 04 Oct 2024
    factor1=sqrt(var1)/sqrt(2);  % scale parameter: b= sigma/sqrt(2);

elseif method1==3 | method1==4 | method1==5

    % Revision SW-1: the synthetic datasets of this ranked model are drawn
    % with ITS OWN fitted dispersion. Before, factor1 kept the value the
    % search left in the global, i.e. the alpha of the top-ranked candidate
    % of the search, for every rank.
    factor1=alpha_hat;

end


plot(timevect,totinc,'r')
hold on
plot(timevect,data,'bo')
xlabel('Time (days)');
ylabel('Incidence')


% <=============================================================================================>
% <======================== Vector with parameter estimates ====================================>
% <=============================================================================================>

Ptrue=[r_hat p_hat a_hat K_hat q_hat alpha_hat d_hat];
fprintf('\nRanked model %d, estimates [r p a K q alpha d]: %s\n',rank1,mat2str(full(Ptrue),5));
fprintf('Bootstrap: %d refits of ranked model %d ...\n',M,rank1);

% <=============================================================================================>
% <========== Generate parameter uncertainty via parametric bootstrapping ======================>
% <=============================================================================================>

yi=cumsum(totinc);

% Synthetic datasets and the start points of their refits, all drawn here on
% the client: the ranked vector, one jittered copy of it and
% fitSettings.bootstrapStarts Latin-hypercube points (default 1), in the
% scaled search space.
curves=zeros(length(yi),M);
startSets=cell(M,1);
unitMap=spatialWaveUnitMap('make',struct('z0',P,'LB',LB,'UB',UB,'isLog',isLog));
for real=1:M
    curves(:,real)=AddPoissonError(yi,1,dist1,factor1,d_hat);
    startSets{real}=spatialWaveStartPoints(unitMap,fitSettings.bootstrapStarts,1,P);
end

[Phatss,bootstrapDiagnostics]=spatialWaveBootstrapRefits(context,curves,P,LB,UB,isLog, ...
    startSets,fitSettings,fitSettings.useParallel);

% <=============================================================================================>
% <================================ Bootstrap summary =========================================>
% <=============================================================================================>

acceptedBootstrap=sum(bootstrapDiagnostics.success);
atLower=sum(abs(Phatss-repmat(LB,M,1))<=1e-6*max(1,abs(repmat(LB,M,1))) & repmat(UB>LB,M,1),1);
atUpper=sum(abs(Phatss-repmat(UB,M,1))<=1e-6*max(1,abs(repmat(UB,M,1))) & repmat(UB>LB,M,1),1);
bootstrapDiagnostics.searchLB=LB;
bootstrapDiagnostics.searchUB=UB;
bootstrapDiagnostics.searchIsLog=isLog;
bootstrapDiagnostics.refitsAtLowerBound=atLower;
bootstrapDiagnostics.refitsAtUpperBound=atUpper;
bootstrapDiagnostics.parameterNames={'r','p','a','K','q','alpha','d'};

% <=========================== Parameter identifiability (PII) ====================================>
% Revision SW-6: the identifiability report of this ranked model
% (spatialWaveIdentifiability.m), from the refits just made: PII, share of
% the search range, share of refits at each bound, Monte Carlo error,
% verdict per parameter and the pairs of parameters that trade off.
bootstrapDiagnostics.identifiability=spatialWaveIdentifiability(Phatss,LB,UB,P,fitSettings);
spatialWavePrintPII(bootstrapDiagnostics.identifiability, ...
    sprintf('Parameter identifiability (PII) of ranked model %d',rank1));
piiReport=bootstrapDiagnostics.identifiability;
notIdentified=piiReport.name(strcmp(piiReport.verdict,'not identified'))';
bootstrapDiagnostics.notIdentified=notIdentified;
if ~isempty(notIdentified)
    warning('SpatialWave:NotIdentified', ...
        ['Ranked model %d: %s not identified by these data. The 95%% interval covers %.0f%% or more of what the ', ...
         'search bounds allow, so it reflects the bounds, not the data; do not interpret the estimate. ', ...
         'See the table printed above.'], ...
        rank1,strjoin(notIdentified,', '),100*piiReport.notIdentifiedRatio);
end

fprintf('\nBootstrap summary for ranked model %d:\n', rank1);
fprintf('Requested bootstrap datasets: %d\n', M);
fprintf('Accepted bootstrap datasets: %d\n', acceptedBootstrap);
fprintf('Unresolved bootstrap datasets: %d\n', M - acceptedBootstrap);
fprintf('Refits with fewer than %d active sub-epidemics: %d\n', nActiveFit, sum(bootstrapDiagnostics.nActive<nActiveFit));
fprintf('Refits with more than %d active sub-epidemics: %d\n', nActiveFit, sum(bootstrapDiagnostics.nActive>nActiveFit));
for k=1:7
    if atLower(k)>0
        fprintf('   %s at its lower bound (%g) in %d of %d refits\n',bootstrapDiagnostics.parameterNames{k},LB(k),atLower(k),M);
    end
    if atUpper(k)>0
        fprintf('   %s at its upper bound (%g) in %d of %d refits\n',bootstrapDiagnostics.parameterNames{k},UB(k),atUpper(k),M);
    end
end
fprintf('Refits in which the polish ran: %d\n\n', sum(bootstrapDiagnostics.polishRun));

% More than 2.5% of the refits at a bound: the 95% percentile interval of
% that parameter ends at the bound, not where the data put it.
for k=[1 4]
    names1={'rLB','rUB';'','';'','';'KLB','KUB'};
    if atLower(k)>0.025*M
        warning('SpatialWave:BootstrapAtBound', ...
            ['Ranked model %d: %d of %d bootstrap refits have %s at its lower limit (%g), so its 95%% interval ', ...
             'ends at the limit. Set %s in options_fit_settings_sw.m to move it.'], ...
            rank1,atLower(k),M,bootstrapDiagnostics.parameterNames{k},LB(k),names1{k,1});
    end
    if atUpper(k)>0.025*M
        warning('SpatialWave:BootstrapAtBound', ...
            ['Ranked model %d: %d of %d bootstrap refits have %s at its upper limit (%g), so its 95%% interval ', ...
             'ends at the limit. Set %s in options_fit_settings_sw.m to move it.'], ...
            rank1,atUpper(k),M,bootstrapDiagnostics.parameterNames{k},UB(k),names1{k,2});
    end
end

% Public point estimates correspond to bestfit and AICc_best.
P0=fitResult.parameters;
d=P0(end);

% Output: the number of ACTIVE sub-epidemics of the accepted fit, as in the
% ranking. fitResult.patchCap holds the number the model was fitted with.
npatches=fitResult.modelSpec(1);
end
