function out = spatialWaveReproductionNumber(S,forecastingperiod,rt)
%SPATIALWAVEREPRODUCTIONNUMBER Rt of one ranked model, with uncertainty.
%   out = spatialWaveReproductionNumber(S,forecastingperiod,rt)
%
%   S                  the contents of the file of a ranked model
%                      (S = load('./output/modifiedLogisticPatch-original-...-rank-k.mat')).
%                      Used: Phatss, P0, data1, patchCap, onset_fixed,
%                      onset_thr, flag1, typedecline1, and fitSettings if present.
%   forecastingperiod  time units beyond the calibration period (0 = none)
%   rt                 struct with the fields type (1 gamma, 2 exponential,
%                      3 delta), mean, variance and minCoverage; [] reads
%                      options_Rt_sw.m.
%
%   HOW. The curve of every bootstrap refit is simulated over the calibration
%   and forecast periods exactly as plotForecast_SW_subepidemicFramework.m
%   simulates it (the sub-epidemics the model was fitted with, and further
%   ones only after the calibration period), and its Rt is computed by the
%   renewal equation (spatialWaveRt.m). The spread over the refits is the
%   uncertainty of Rt that comes from the parameter estimates. The Rt of the
%   ranked fit itself (S.P0) is returned as the estimate.
%
%   Rt is computed from the model curve, not from the data, and without
%   observation noise. If the fit used the day-of-week adjustment (revision
%   SW-5), the curve WITHOUT the day-of-week factors is used: a weekly
%   reporting cycle is not a weekly cycle of transmission.
%
%   OUTPUT, a struct:
%     time        whole times 1, 2, ..., T (calibration and forecast), T-by-1
%     calibration true for the times of the calibration period
%     estimate    Rt of the ranked fit at those times
%     median, LB, UB   median and 95% interval over the refits
%     aboveOne    share of the refits with Rt above 1 at each time
%     nRefits     refits used; nSkipped: refits whose curve could not be
%                 simulated (left out; see firstSkipMessage)
%     firstTime   first time with an Rt value (earlier ones are blank
%                 because the curve does not yet cover the generation
%                 interval); NaN if there is none
%     atCalibrationEnd   [time estimate median LB UB aboveOne] at the last
%                 calibration time
%     atEnd       the same at the last time
%     ts, Rss     the fine grid (step 0.01) and Rt of every refit on it:
%                 Rss(k,b) is the value of refit b at time ts(k+1)
%     rt          the generation-interval settings used
%   Times before firstTime hold NaN.
%
%   READING IT. Rt here is what the fitted curve implies under the stated
%   generation interval. It inherits every limit of the fit: where the model
%   forecasts badly, its Rt beyond the calibration period is wrong too. The
%   interval reflects parameter uncertainty only, given the model.
%
%   Revision SW-7 (after plotReproductionNumber.m of the SubEpiPredict toolbox).
%
%   See also spatialWaveRt, plotReproductionNumber_SW, options_Rt_sw

if nargin<2 || isempty(forecastingperiod)
    forecastingperiod=0;
end
if ~(isnumeric(forecastingperiod) && isscalar(forecastingperiod) && forecastingperiod>=0 && forecastingperiod==fix(forecastingperiod))
    error('SpatialWave:InvalidRtSettings','forecastingperiod must be a whole number, 0 or more.');
end
if nargin<3 || isempty(rt)
    [typeGI,meanGI,varGI,minCoverage]=options_Rt_sw;
    rt=struct('type',typeGI,'mean',meanGI,'variance',varGI,'minCoverage',minCoverage);
end
needed={'Phatss','P0','data1','patchCap','onset_fixed','onset_thr','flag1','typedecline1'};
for k=1:numel(needed)
    if ~isfield(S,needed{k})
        error('SpatialWave:InvalidRtInput', ...
            'The fit does not hold %s. Rt needs the file of a ranked model saved by revision SW-1 or later.',needed{k});
    end
end
data1=S.data1;
nCal=size(data1,1);
if nCal<2
    error('SpatialWave:InvalidRtInput','The calibration period has fewer than 2 observations.');
end
timevect2=(0:nCal-1+forecastingperiod)';

% the model without the day-of-week factors (see the header)
fitSettings=[];
if isfield(S,'fitSettings')
    fitSettings=S.fitSettings;
end
weekdayCurve=spatialWaveWeekdayCurve(spatialWaveWeekdayFactors(data1(:,2),fitSettings),1);
if isempty(weekdayCurve)
    I0=data1(1,2);
else
    I0=data1(1,2)/weekdayCurve(1);
end
waveExtension=struct('extraPatches',spatialWaveForecastExtraPatches(),'unlockTime',data1(end,1));

M=size(S.Phatss,1);
Rss=[];
ts=[];
nSkipped=0;
firstSkipMessage='';
for b=0:M
    if b==0
        z=S.P0;          % the ranked fit
    else
        z=S.Phatss(b,:);
    end
    try
        [~,~,totinc]=simulateSpatialWave(timevect2,I0,S.patchCap,S.onset_fixed,S.onset_thr, ...
            S.flag1,S.typedecline1,z(1),z(2),z(3),z(4),z(5),[],waveExtension);
    catch simulationError
        if ~strncmp(simulationError.identifier,'SpatialWave:',12) || b==0
            rethrow(simulationError);
        end
        nSkipped=nSkipped+1;
        if isempty(firstSkipMessage)
            firstSkipMessage=simulationError.message;
        end
        continue
    end
    [Rs,ts]=spatialWaveRt(timevect2,totinc,rt.type,rt.mean,rt.variance,rt.minCoverage);
    if b==0
        Rfit=Rs;
    else
        Rss=[Rss Rs]; %#ok<AGROW>
    end
end

% one row per whole time 1, 2, ..., T: the value at whole time j is entry 100*j
T=numel(timevect2)-1;
rowsRt=round((1:T)'*100);
out=struct();
out.time=(1:T)';
out.calibration=out.time<=nCal-1;
out.estimate=Rfit(rowsRt);
if isempty(Rss)
    out.median=NaN(T,1); out.LB=NaN(T,1); out.UB=NaN(T,1); out.aboveOne=NaN(T,1);
else
    whole=Rss(rowsRt,:);
    out.median=localQuantile(whole,0.5);
    out.LB=localQuantile(whole,0.025);
    out.UB=localQuantile(whole,0.975);
    out.aboveOne=sum(whole>1,2)./max(sum(isfinite(whole),2),1);
    out.aboveOne(~any(isfinite(whole),2))=NaN;
end
out.nRefits=size(Rss,2);
out.nSkipped=nSkipped;
out.firstSkipMessage=firstSkipMessage;
first=find(isfinite(out.estimate),1,'first');
if isempty(first)
    out.firstTime=NaN;
else
    out.firstTime=out.time(first);
end
rowCal=nCal-1;
out.atCalibrationEnd=[out.time(rowCal) out.estimate(rowCal) out.median(rowCal) out.LB(rowCal) out.UB(rowCal) out.aboveOne(rowCal)];
out.atEnd=[out.time(T) out.estimate(T) out.median(T) out.LB(T) out.UB(T) out.aboveOne(T)];
out.ts=ts;
out.Rss=Rss;
out.rt=rt;
end

% -----------------------------------------------------------------------------
function q = localQuantile(A,p)
% Row-wise quantile that leaves a row of NaN as NaN.
q=NaN(size(A,1),1);
for k=1:size(A,1)
    x=A(k,isfinite(A(k,:)));
    if ~isempty(x)
        q(k)=quantile(x,p);
    end
end
end
