function report=spatialWaveSelfCheck()
%SPATIALWAVESELFCHECK Quick checks of the simulation layer (SW-1) and the fit layer (SW-2).
%   report = spatialWaveSelfCheck() runs in MATLAB (and in GNU Octave) with
%   the toolbox folder on the path, stops with an error at the first check
%   that fails and otherwise prints one line per check. No files and no
%   figures are involved; one small model is fitted (check 13). It takes
%   some seconds.
%
%   Checks 2 and 3 are the ones revision SW-1 could not run in MATLAB: they
%   ask that an activation is located precisely on the UNIT time grid the
%   toolbox fits on. Octave's ode15s locates events between output times
%   only, so the tolerance is loose there; in MATLAB it is 1e-6.
%
%   runtests('tests') runs this function through test_spatialWaveSelfCheck.m.

isOctave=exist('OCTAVE_VERSION','builtin')~=0;
if isOctave
    eventTol=1e-2;   % time units, unit grid
    stateTol=5e-3;   % relative, unit grid against fine grid
else
    eventTol=1e-6;
    stateTol=1e-6;
end
report={};

% 1. one logistic sub-epidemic against its closed form
r=0.25; K=1000; I0=5; t=(0:90)';
[~,x]=simulateSpatialWave(t,I0,1,0,0,3,1,r,1,1,K,0.3);
exact=K./(1+(K/I0-1)*exp(-r*t));
check(max(abs(x-exact)./exact)<1e-6,'logistic sub-epidemic equals its closed form');

% 2. activation time on the unit grid against its closed form
thr=400;
[~,~,~,~,act]=simulateSpatialWave((0:60)',I0,3,0,thr,3,1,r,1,1,K,0.3);
expected=log((K/I0-1)/(K/thr-1))/r;
check(abs(act.time(2)-expected)<eventTol, ...
    sprintf('activation time on the unit grid (error %.2g, limit %.0g)',abs(act.time(2)-expected),eventTol));

% 3. the unit grid and a grid 100 times finer give the same states
args={5,4,0,400,1,1,0.3,0.9,1,1000,0.3};
[~,xUnit,~,~,actUnit]=simulateSpatialWave((0:90)',args{:});
[~,xFine,~,~,actFine]=simulateSpatialWave((0:0.01:90)',args{:});
difference=max(max(abs(xUnit-xFine(1:100:end,:))))/max(xUnit(:));
check(difference<stateTol,sprintf('unit grid equals fine grid (difference %.2g, limit %.0g)',difference,stateTol));
check(actUnit.nActive==3 && actFine.nActive==3,'three sub-epidemics active after 90 time units');
check(abs(actFine.time(2)-24.364377)<1e-4 && abs(actFine.time(3)-56.068874)<1e-4, ...
    'activation times 24.3644 and 56.0689 (values of the Python reference)');

% 4. incidence conventions
[~,x,totinc,patchIncidence]=simulateSpatialWave((0:90)',args{:});
check(abs(totinc(1)-5)<1e-12,'first incidence is the first observation');
check(max(abs(sum(patchIncidence,2)-totinc))<1e-9,'sub-epidemic incidences add up to the total');
check(all(x(:,4)==1),'a sub-epidemic that has not started keeps its seed');

% 5. a forecast curve equals the fitted curve over the calibration period
extension=struct('extraPatches',10,'unlockTime',60);
[~,xFit]=simulateSpatialWave((0:60)',5,2,0,200,1,1,0.3,0.9,1,1000,0.1);
[~,xForecast,~,~,actForecast]=simulateSpatialWave((0:150)',5,2,0,200,1,1,0.3,0.9,1,1000,0.1,[],extension);
check(max(max(abs(xForecast(1:61,1:2)-xFit)))/max(xFit(:))<max(1e-9,stateTol*1e-3), ...
    'forecast equals fit over the calibration period');
check(size(xForecast,2)==12 && all(all(xForecast(1:61,3:end)==1)),'extra sub-epidemics dormant until the unlock time');
check(actForecast.immediateAtUnlock && actForecast.time(3)==60,'a binding limit is reported at the unlock time');

% 6. synchronous start
[~,x,totinc]=simulateSpatialWave((0:30)',6,3,1,0,1,2,0.3,0.9,1,1000,0.7);
check(all(abs(x(1,:)-2)<1e-12) && abs(totinc(1)-6)<1e-12,'synchronous start from I0/n each');

% 7. a threshold above K never starts a second sub-epidemic
[~,~,~,~,act]=simulateSpatialWave((0:90)',5,4,0,1500,1,1,0.3,0.9,1,1000,0.3);
check(act.nActive==1,'threshold above K: one sub-epidemic');
[~,~,~,~,act]=simulateSpatialWave((0:20)',5,3,0,20,0,1,0.3,0.6,1,1,0);
check(act.nActive==1,'GGM with a threshold above 1: one sub-epidemic');

% 8. the reported wave is the wave the dynamics produce
[Ks,Ktot,n]=getsubepidemicsizesFinite(1000,400,0.3,1,14);
[~,x,~,~,act]=simulateSpatialWave((0:5:1500)',5,14,0,400,3,1,0.6,1,1,1000,0.3);
check(act.nActive==n,sprintf('wave of %d sub-epidemics, as getsubepidemicsizesFinite says',n));
check(abs(sum(x(end,1:n))-Ktot)/Ktot<1e-6 && numel(Ks)==n,'wave size equals the sum of the sizes');

% 9. negative-binomial objective against the original loops
y=[0 3 7 12 40 5]'; mu=[0.4 2.5 9 10 33 6.5]';
for method=3:5
    alpha=0.7; d=1.3; sum1=0;
    for i=1:numel(y)
        for j=0:(y(i)-1)
            switch method
                case 3, sum1=sum1+log(j+(1/alpha)*mu(i));
                case 4, sum1=sum1+log(j+(1/alpha));
                case 5, sum1=sum1+log(j+(1/alpha)*mu(i)^(2-d));
            end
        end
        switch method
            case 3, sum1=sum1+y(i)*log(alpha)-(y(i)+(1/alpha)*mu(i))*log(1+alpha);
            case 4, sum1=sum1+y(i)*log(alpha*mu(i))-(y(i)+(1/alpha))*log(1+alpha*mu(i));
            case 5, sum1=sum1+y(i)*log(alpha*mu(i)^(d-1))-(y(i)+(1/alpha)*mu(i)^(2-d))*log(1+alpha*mu(i)^(d-1));
        end
    end
    value=getSpatialWaveObjective(y,mu,method,alpha,d);
    check(abs(value+sum1)<1e-9*max(1,abs(sum1)),sprintf('negative-binomial objective, method %d',method));
end

% ---- revision SW-2: the fit layer -------------------------------------------

% 10. unit coordinates: every free parameter maps to [0,1] and back
settings=spatialWaveFitSettings();
settings.useParallel=false; settings.polish=false; settings.fitMaxAttempts=1;
t=(0:40)'; truth=[0.3 0.9 1 1000 0.3 0 0];
[~,~,y]=simulateSpatialWave(t,5,2,1,0,1,1,truth(1),truth(2),truth(3),truth(4),truth(5));
[LB,UB,isLog,z0]=spatialWaveSearchSpace(y,1,0,settings);
map=spatialWaveUnitMap('make',struct('z0',z0,'LB',LB,'UB',UB,'isLog',isLog));
u=spatialWaveUnitMap('toU',truth,map);
back=spatialWaveUnitMap('toZ',u,map);
check(map.nFree==4 && all(u>=0 & u<=1) && max(abs(back-truth)./max(1,abs(truth)))<1e-10, ...
    'unit coordinates: four free parameters, there and back');
check(isLog(1) && isLog(4) && ~isLog(2) && UB(1)>=10*max(y), ...
    'r and K on a log scale; upper bound of r follows the counts');

% 11. the objective with an explicit context is zero at the true parameters
context=struct('timevect',t,'ydata',y,'I0',5,'npatches',2,'onset_fixed',1,'onset_thr',0, ...
    'flag1',1,'typedecline1',1,'method1',0,'odeOptions',spatialWaveOdeOptions());
[value,detail]=spatialWaveObjective(truth,context);
check(detail.success && value<1e-12 && detail.nActive==2,'context objective is zero at the true parameters');

% 12. thresholds: the low end of the grid is filled, on a log scale
counts=[3 3 5 11 18 37 43 70 61 53 39 34 40 46 45 39 43 33 20 20]';
coarse=spatialWaveThresholdGrid(counts,1,3,0);
fine=spatialWaveThresholdGrid(counts,1,3,6);
added=setdiff(fine,coarse);
check(numel(coarse)==20 && numel(added)==6 && all(added>coarse(2) & added<coarse(3)) && ...
    max(abs(diff(log(added))-mean(diff(log(added)))))<1e-9,'six thresholds added below the first step of the grid');

% ---- revision SW-3: flat baseline, rolling summary, fit signature ------------
base=spatialWaveBaselineForecast((1:20)',2,1);
check(base.level==20 && isequal(base.nErrors,[19;18]) && abs(median(base.samples(1,:))-20)<1e-12 && ...
    abs(min(base.samples(1,:))-19)<1e-12 && abs(max(base.samples(1,:))-21)<1e-12 && ...
    abs(min(base.samples(2,:))-18)<1e-12 && abs(max(base.samples(2,:))-22)<1e-12, ...
    'flat baseline: last value, with the past errors of the same rule in both directions');
base=spatialWaveBaselineForecast([zeros(10,1)+4;6;8],3,7);
check(abs(base.level-(5*4+6+8)/7)<1e-12 && all(base.samples(:)>=0) && all(all(diff(base.samples,1,2)>=0)), ...
    'flat baseline: mean of the last 7 values, samples in order and not below 0');
% two dates, three steps; the model has half the WIS of the baseline at every step
rollScores=zeros(0,8);
for rollDate=[1 2]
    for rollStep=1:3
        baseWIS=10*rollDate;                         % the same at every step, so also over the first h steps
        rollScores=[rollScores; rollDate 0 rollStep 2*baseWIS 0 100 baseWIS baseWIS; ...
                                rollDate 1 rollStep baseWIS 0 100 baseWIS/2 baseWIS]; %#ok<AGROW>
    end
end
[rollSummary,rollColumns]=spatialWaveRollingSummary(rollScores,3);
rowModel=rollSummary(rollSummary(:,1)==1 & rollSummary(:,2)==3,:);
check(numel(rollColumns)==13 && size(rollSummary,1)==6 && rowModel(3)==2 && abs(rowModel(6)-7.5)<1e-12 && ...
    abs(rowModel(7)-0.5)<1e-12 && abs(rowModel(8)-0.5)<1e-12 && rowModel(9)==1 && abs(rowModel(10)-0.5)<1e-12 && ...
    abs(rowModel(11)-7.5)<1e-12 && abs(rowModel(12)-0.5)<1e-12, ...
    'rolling summary: a model with half the WIS of the baseline gets ratio 0.5 and share 1');
savedFit=struct('fitRevision',spatialWaveFitRevision(),'fitSettings',spatialWaveFitSettings());
sameA=spatialWaveFitSignature(savedFit);
savedFit.fitSettings.useParallel=~savedFit.fitSettings.useParallel;
sameB=spatialWaveFitSignature(savedFit);
savedFit.fitSettings.thresholdLowPoints=savedFit.fitSettings.thresholdLowPoints+1;
[sameC,whyC]=spatialWaveFitSignature(savedFit);
sameD=spatialWaveFitSignature(struct('M',300));
check(sameA && sameB && ~sameC && numel(whyC)==1 && ~sameD, ...
    'fit signature: same settings reuse a fit, another setting or an unrecorded fit does not');

% ---- revision SW-5: day-of-week factors ----------------------------------------
cycle=[1.3;1.2;1.1;1.0;0.9;0.8;0.7];
weekdayOn=struct('weekdayAdjustment',true);
wd=spatialWaveWeekdayFactors(100*repmat(cycle,4,1),weekdayOn);
wdOff=spatialWaveWeekdayFactors(100*repmat(cycle,4,1),struct('weekdayAdjustment',false));
wdCurve=spatialWaveWeekdayCurve(wd,16);
check(wd.applied && max(abs(wd.byPosition-cycle))<1e-12 && isequal(wd.nRatios,[3;3;3;4;3;3;3]) && ...
    isequal(wdCurve,cycle([1:7 1:7 1 2])) && ~wdOff.applied && isempty(spatialWaveWeekdayCurve(wdOff,16)), ...
    'day-of-week factors of an exact weekly cycle; none when switched off');
holiday=100*repmat(cycle,6,1); holiday(18)=5;                 % one day far too low
wdHoliday=spatialWaveWeekdayFactors(holiday,weekdayOn);
check(max(abs(wdHoliday.byPosition-cycle))<0.02,'day-of-week factors: one unusual day does not move them');
% the curve of check 11 with a weekly cycle on it: the objective with the factors is zero at the truth
weekCurve=spatialWaveWeekdayCurve(wd,numel(y));
contextWeek=context; contextWeek.ydata=y(:).*weekCurve; contextWeek.weekdayFactor=weekCurve;
[valueWeek,detailWeek]=spatialWaveObjective(truth,contextWeek);
contextPlain=contextWeek; contextPlain.weekdayFactor=[];
valuePlain=spatialWaveObjective(truth,contextPlain);
check(detailWeek.success && valueWeek<1e-12*sum(y.^2) && valuePlain>1e-3*sum(y.^2), ...
    'objective with day-of-week factors is zero at the true parameters');
savedFit=struct('fitRevision',spatialWaveFitRevision(),'fitSettings',rmfield(spatialWaveFitSettings(),'weekdayAdjustment'));
nowOn=spatialWaveFitSettings();
check(spatialWaveFitSignature(savedFit)==~nowOn.weekdayAdjustment, ...
    'a fit saved before the day-of-week setting existed is usable exactly while the setting is off');

% ---- revision SW-6: identifiability report (PII) --------------------------------
[piiA,fracA,scaleA,boxA,ratioA]=spatialWavePII([10 0.2],[1000 0.5],[1 0],[1e4 1]);
check(abs(piiA(1)-2)<1e-12 && abs(boxA(1)-4)<1e-12 && abs(ratioA(1)-0.5)<1e-12 && isequal(scaleA,[1 0]) && ...
    abs(fracA(2)-0.3)<1e-12 && abs(ratioA(2)-0.3)<1e-12 && isnan(boxA(2)), ...
    'PII: 10 to 1000 is 2 decades, half of what bounds 1 to 1e4 allow; a lower bound of 0 uses the range');
[seLow,seHigh,sePII]=spatialWaveQuantileMCSE((1:50)',20);
check(abs(seLow-1.4078609194167768)<1e-9 && abs(seHigh-1.1323589722809915)<1e-9 && abs(sePII-0.25366590555870855)<1e-9, ...
    'Monte Carlo error of the interval (values of the Python reference)');
ridge=linspace(-1,1,200)'; wobble=sin(7*(1:200)');
refits=[0.3*exp(0.02*wobble) 0.9+0.005*cos(11*(1:200)') ones(200,1) 1e4*10.^ridge 1.5+1.4*ridge ones(200,2)];
reportSettings=struct('PIIscale',[],'notIdentifiedRatio',0.5,'weaklyIdentifiedPII',1,'tradeoffCorrelation',0.9);
piiReport=spatialWaveIdentifiability(refits,[1e-3 0 1 20 0 1 1],[500 1 1 1e8 5 1 1],median(refits,1),reportSettings);
check(isequal(piiReport.name',{'r','p','a','K','q'}) && ...
    isequal(piiReport.verdict',{'identified','identified','fixed','weakly identified','not identified'}) && ...
    piiReport.nFree==4 && piiReport.nIdentified==2 && abs(piiReport.PII(4)-1.9)<0.02 && abs(piiReport.PII_ratio(5)-0.532)<0.005, ...
    'identifiability: K over 2 decades is weakly identified, q over half its range is not');
check(numel(piiReport.tradeoffs)==1 && strcmp(piiReport.tradeoffs{1},'K-q (+1.00)') && abs(piiReport.correlation(4,5)-1)<1e-12 && ...
    abs(piiReport.correlation(1,4))<0.5 && all(isnan(piiReport.correlation(3,:))), ...
    'identifiability: K and q, which move together over the refits, are reported as a pair');
check(isempty(spatialWaveSavedIdentifiability(refits,median(refits,1),struct('nActive',1))) && ...
    ~isempty(spatialWaveSavedIdentifiability(refits,median(refits,1),struct('searchLB',[1e-3 0 1 20 0 1 1],'searchUB',[500 1 1 1e8 5 1 1]))), ...
    'identifiability of a saved fit: from its refits and search bounds, none without the bounds');

% ---- revision SW-7: effective reproduction number --------------------------------
rtTimes=(0:60)';
[rtFlat,~,rtCover]=spatialWaveRt(rtTimes,100+zeros(61,1),1,5,4);
check(isnan(rtFlat(100)) && rtCover(100)<0.99 && abs(rtFlat(end)-1)<0.011 && all(isnan(rtFlat(rtCover<0.99-1e-9))), ...
    'Rt of a constant incidence is 1, and blank until the curve covers the generation interval');
rtUp=spatialWaveRt(rtTimes,50*exp(0.08*rtTimes),1,5,4);
rtDown=spatialWaveRt(rtTimes,50*exp(-0.05*rtTimes),3,6,[]);
check(abs(rtUp(end)-(1+0.08*0.8)^(5/0.8))<5e-3*(1+0.08*0.8)^(5/0.8) && abs(rtDown(end)-exp(-0.05*6))<2e-3, ...
    'Rt of exponential growth equals its closed form (gamma and fixed generation interval)');
rtTruth=[0.6 0.9 1 3000 0.4 1 1];
[~,~,rtCurve]=simulateSpatialWave((0:39)',4,1,0,1e6,1,1,rtTruth(1),rtTruth(2),rtTruth(3),rtTruth(4),rtTruth(5));
rtRefits=repmat(rtTruth,40,1); rtRefits(:,1)=rtTruth(1)*(1+0.05*sin((1:40)')); rtRefits(:,4)=rtTruth(4)*(1+0.1*cos((1:40)'));
rtFit=struct('Phatss',rtRefits,'P0',rtTruth,'data1',[(0:39)' rtCurve(:)],'patchCap',1,'onset_fixed',0,'onset_thr',1e6,'flag1',1,'typedecline1',1);
rtOut=spatialWaveReproductionNumber(rtFit,10,struct('type',1,'mean',3,'variance',2,'minCoverage',0.99));
rtShown=isfinite(rtOut.median);
check(numel(rtOut.time)==49 && sum(rtOut.calibration)==39 && rtOut.nRefits==40 && rtOut.nSkipped==0 && rtOut.firstTime==8 && ...
    all(rtOut.LB(rtShown)<=rtOut.median(rtShown)+1e-12) && all(rtOut.median(rtShown)<=rtOut.UB(rtShown)+1e-12) && ...
    rtOut.estimate(8)>1 && rtOut.atEnd(2)<1 && rtOut.atCalibrationEnd(1)==39 && all(diff(rtOut.estimate(rtShown))<0), ...
    'Rt of a fitted wave: above 1 while it grows, below 1 after it has turned, interval around the median');

% ---- revision SW-8: quantile scores and intervals that allow for model error ------
[~,~,qLevels]=spatialWaveQuantileScores(zeros(1,23),0);
qAlphas=[0.02 0.05 0.1:0.1:0.9];
[wisIn,insideIn]=spatialWaveQuantileScores(100*qLevels,50);          % quantiles of a uniform forecast on 0 to 100
[wisOut,insideOut]=spatialWaveQuantileScores(100*qLevels,200);
wisInHand=sum((qAlphas/2).*(100*(1-qAlphas)))/11.5;
wisOutHand=(0.5*150+sum((qAlphas/2).*(100*(1-qAlphas)+(2./qAlphas).*(200-100*(1-qAlphas/2)))))/11.5;
check(numel(qLevels)==23 && insideIn==1 && insideOut==0 && abs(wisIn-wisInHand)<1e-10 && abs(wisOut-wisOutHand)<1e-10, ...
    'WIS and coverage from quantiles: an observation inside and one outside the forecast');
% one model, 12 daily forecast dates, one step: intervals of 100 +/- 10, observations 30 away
qNormal=[-2.3263 -1.96 -1.6449 -1.2816 -1.0364 -0.8416 -0.6745 -0.5244 -0.3853 -0.2533 -0.1257 0 ...
    0.1257 0.2533 0.3853 0.5244 0.6745 0.8416 1.0364 1.2816 1.6449 1.96 2.3263];
errFc=zeros(12,27);
for errDate=1:12
    errFc(errDate,:)=[errDate 1 1 100+30*(-1)^errDate 100+(10/1.96)*qNormal];
end
errFc(11:12,4)=[80;120];   % the last two observations 20 away: clear of the widened limits of 70 and 130
[errWide,errScales]=spatialWaveModelErrorIntervals(errFc,1,struct('modelErrorWindow',3,'modelErrorMinErrors',10,'modelErrorMaxScale',20));
check(isequal(errScales(:,5)',0:11) && all(errScales(1:10,4)==1) && all(abs(errScales(11:12,4)-3)<1e-9) && ...
    isequal(errWide(1:10,:),errFc(1:10,:)) && abs(errWide(12,4+2)-70)<1e-6 && abs(errWide(12,4+22)-130)<1e-6 && errWide(12,4+12)==100, ...
    'model-error intervals: errors 3 times the half-width widen the interval 3 times, from 10 past errors on');
errBefore=spatialWaveForecastScores(errFc,[errFc(:,1:3) zeros(12,5)]);
errAfter=spatialWaveForecastScores(errWide,[errFc(:,1:3) zeros(12,5)]);
check(all(errBefore(:,6)==0) && all(errAfter(11:12,6)==100) && all(errAfter(11:12,7)<errBefore(11:12,7)) && ...
    isequal(errAfter(1:10,6:7),errBefore(1:10,6:7)), ...
    'model-error intervals: coverage rises from 0 to 100% and the WIS falls where the interval was widened');

% 13. the fit layer finds a fit at least as good as the true parameters.
%     In MATLAB this runs MultiStart and fmincon with the default settings;
%     it is the first check that does. The series is the model curve with a
%     fixed +/-15% disturbance, like real data. The criterion is the
%     objective, not the parameters: an optimizer that has converged does at
%     least as well as the truth on the series it was given. (The first
%     version of this check asked for the parameters of a curve WITHOUT
%     noise to 1%, with the polish switched off. That is a poor test: the
%     scaled objective, log(SSE), has no lower bound on such a series, and
%     how close the parameters come also depends on how well they are
%     identified. tests/spatialWaveFitDiagnostic.m shows both.)
%     Octave has neither MultiStart nor fmincon: the check is skipped there
%     unless stand-ins are on the path, with few start points and a loose limit.
k=(1:numel(t))';
noisy=y.*(1+0.15*sin(1.7*k)); noisy(1)=y(1);
fitContext=context; fitContext.ydata=noisy;
fitSettings=spatialWaveFitSettings(); fitSettings.useParallel=false;
if isOctave
    fitSettings.polish=false; fitSettings.fitMaxAttempts=1; nRandomStarts=1; objectiveLimit=1.10;
else
    nRandomStarts=5; objectiveLimit=1+1e-6;
end
[LBn,UBn,isLogn,z0n]=spatialWaveSearchSpace(noisy,1,0,fitSettings);
sseTruth=spatialWaveObjective(truth,fitContext);
canFit=true;
try
    [P,fval,fitInfo]=fitSpatialWaveModel(fitContext,z0n,LBn,UBn,isLogn,[],nRandomStarts,fitSettings);
catch fitError
    if isOctave
        canFit=false;
        fprintf('  --  fit check skipped: %s\n',strtrim(fitError.message));
    else
        rethrow(fitError);
    end
end
if canFit
    check(fitInfo.success && fitInfo.nActive==2, ...
        sprintf('fit accepted (%d local runs, %d converged)',fitInfo.nLocalRuns,fitInfo.nConvergedRuns));
    check(fval<=objectiveLimit*sseTruth, ...
        sprintf('fit at least as good as the true parameters (SSE %.5g against %.5g, limit %.7g times)', ...
        fval,sseTruth,objectiveLimit));
    check(P(3)==1 && P(6)==0 && P(7)==0,'fixed parameters keep their value');
    fprintf('      fitted [r p K q] = %s, true %s\n',mat2str(P([1 2 4 5]),4),mat2str(truth([1 2 4 5]),4));
end

fprintf('spatialWaveSelfCheck: %d checks passed.\n',numel(report));

    function check(condition,label)
        if ~condition
            error('SpatialWave:SelfCheckFailed','Check failed: %s',label);
        end
        report{end+1}=label;
        fprintf('  ok  %s\n',label);
    end
end
